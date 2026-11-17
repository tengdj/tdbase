#include "api.h"
#include <chrono>
#include "include/farm.h"
#include <filesystem>
#include <cstdlib>
#include <memory>
#include <stdexcept>

namespace tdbase::farm2d {
namespace {

// Bound every read before handing WKT to the original geometry parser.
void validate_wkt(const std::string &line, std::size_t row)
{
    std::size_t at = line.find_first_not_of(" \t\r");
    const auto fail = [&] { throw std::invalid_argument("Invalid POLYGON WKT at line " + std::to_string(row)); };
    if(at == std::string::npos || line.compare(at, 7, "POLYGON") != 0) fail();
    at += 7;
    const auto space = [&] { while(at < line.size() && std::isspace(static_cast<unsigned char>(line[at]))) ++at; };
    const auto token = [&](char value) { space(); if(at >= line.size() || line[at++] != value) fail(); };
    const auto number = [&] {
        space();
        if(at >= line.size()) fail();
        char *end = nullptr;
        const double value = std::strtod(line.c_str() + at, &end);
        if(end == line.c_str() + at || !std::isfinite(value)) fail();
        at = static_cast<std::size_t>(end - line.c_str());
        return value;
    };
    token('(');
    std::size_t rings = 0;
    for(;;){
        if(++rings > 1)
            throw std::invalid_argument("FARM query kernels require POLYGON without interior rings (line " + std::to_string(row) + ")");
        token('(');
        std::size_t vertices = 0;
        double first_x = 0, first_y = 0, last_x = 0, last_y = 0;
        double min_x = HUGE_VAL, max_x = -HUGE_VAL, min_y = HUGE_VAL, max_y = -HUGE_VAL;
        for(;;){
            last_x = number(); last_y = number();
            if(!vertices){ first_x = last_x; first_y = last_y; }
            ++vertices;
            min_x = std::min(min_x, last_x); max_x = std::max(max_x, last_x);
            min_y = std::min(min_y, last_y); max_y = std::max(max_y, last_y);
            space();
            if(at >= line.size()) fail();
            if(line[at] != ',') break;
            ++at;
        }
        token(')');
        if(vertices < 4 || first_x != last_x || first_y != last_y || min_x >= max_x || min_y >= max_y) fail();
        space();
        if(at >= line.size()) fail();
        if(line[at] != ',') break;
        ++at;
    }
    token(')'); space();
    if(at != line.size()) fail();
}

// Legacy FARM .dt files have a metadata footer, unlike the 3D tile format.
std::size_t validate_binary(const std::string &path)
{
    std::ifstream in(path, std::ios::binary | std::ios::ate);
    if(!in) throw std::runtime_error("Cannot open 2D dataset: " + path);
    const auto length = in.tellg();
    if(length < static_cast<std::streamoff>(sizeof(std::size_t)))
        throw std::runtime_error("Invalid FARM dataset: " + path);
    std::size_t count = 0;
    in.seekg(-static_cast<std::streamoff>(sizeof(count)), std::ios::end);
    in.read(reinterpret_cast<char *>(&count), sizeof(count));
    const auto bytes = static_cast<std::size_t>(length) - sizeof(count);
    if(count > bytes / sizeof(PolygonMeta))
        throw std::runtime_error("Expected a FARM 2D dataset: " + path);
    const auto payload = bytes - count * sizeof(PolygonMeta);
    in.seekg(static_cast<std::streamoff>(payload));
    std::size_t previous_end = 0;
    for(std::size_t i = 0; i < count; ++i){
        PolygonMeta meta{};
        in.read(reinterpret_cast<char *>(&meta), sizeof(meta));
        if(!in || meta.offset != previous_end || meta.offset > payload ||
           meta.size > payload - meta.offset || meta.num_vertices < 4 ||
           meta.size < sizeof(int) + meta.num_vertices * sizeof(Point))
            throw std::runtime_error("Invalid FARM polygon metadata: " + path);
        previous_end = meta.offset + meta.size;
    }
    if(previous_end != payload)
        throw std::runtime_error("Invalid FARM dataset length: " + path);
    return count;
}

// GPU kernels use compact array positions; public results keep dataset IDs.
std::unique_ptr<Farm> load_source_object(const std::string &path, std::size_t count, std::size_t id)
{
    std::ifstream in(path, std::ios::binary);
    PolygonMeta meta{};
    const auto footer_bytes = sizeof(std::size_t) + (count - id) * sizeof(PolygonMeta);
    in.seekg(-static_cast<std::streamoff>(footer_bytes), std::ios::end);
    if(!in.read(reinterpret_cast<char *>(&meta), sizeof(meta)))
        throw std::runtime_error("Cannot read source object metadata: " + path);
    std::vector<char> buffer(meta.size);
    in.seekg(static_cast<std::streamoff>(meta.offset));
    if(!in.read(buffer.data(), static_cast<std::streamsize>(buffer.size())))
        throw std::runtime_error("Cannot read source object: " + path);
    auto object = std::make_unique<Farm>();
    if(object->decode(buffer.data()) != meta.size || object->get_num_vertices() != meta.num_vertices)
        throw std::runtime_error("Invalid source object: " + path);
    object->getMBB();
    object->id = 0;
    return object;
}

struct CandidateIndex {
    RTree<Farm *, double, 2, double> tree;
    std::int64_t selected_source = -1;
    bool self_join = false;
};

bool collect_candidate(Farm *source, void *arg)
{
    auto *ctx = static_cast<query_context *>(arg);
    auto *global = ctx->global_ctx;
    auto *target = static_cast<Farm *>(ctx->target);
    const auto *index = static_cast<CandidateIndex *>(ctx->target2);
    const auto source_id = index->selected_source >= 0
        ? static_cast<std::size_t>(index->selected_source) : source->id;
    if(index->self_join && source_id == target->id) return true;
    const auto target_id = target->id +
        (global->use_gpu && !global->target_path.empty() ? global->source_objects.size() : 0);
    ctx->object_pairs.emplace_back(source->id, target_id);
    return true;
}

void *collect_candidates(void *arg)
{
    auto *ctx = static_cast<query_context *>(arg);
    auto *global = ctx->global_ctx;
    auto *index = static_cast<CandidateIndex *>(ctx->target2);
    auto &targets = global->target_path.empty() ? global->source_objects : global->target_objects;
    while(ctx->next_batch(10)){
        for(auto i = ctx->index; i < ctx->index_end; ++i){
            Farm *target = targets[i];
            box bounds = global->is_within_query()
                ? target->getMBB()->expand(global->within_distance, true) : *target->getMBB();
            ctx->target = target;
            index->tree.Search(bounds.low, bounds.high, collect_candidate, ctx);
        }
    }
    ctx->merge_object_pairs();
    return nullptr;
}

void *compute_cpu(void *arg)
{
    auto *ctx = static_cast<query_context *>(arg);
    auto *global = ctx->global_ctx;
    auto &targets = global->target_path.empty() ? global->source_objects : global->target_objects;
    while(ctx->next_batch(10)){
        for(auto i = ctx->index; i < ctx->index_end; ++i){
            const auto ids = global->object_pairs[i];
            Farm *source = global->source_objects[ids.first];
            Farm *target = targets[ids.second];
            bool matched = false;
            if(global->query_type == QueryType::within){
                matched = source->within(target, ctx, global->use_approximation);
            }else{
                if(source->get_step_x() < target->get_step_x()) std::swap(source, target);
                if(global->query_type == QueryType::intersection){
                    ctx->target_id = i;
                    source->intersection(target, ctx, global->use_approximation);
                    matched = global->areas[i] > 0;
                }else{
                    matched = source->intersect(target, ctx, global->use_approximation);
                }
            }
            ctx->found += matched;
            if(global->collect_results) global->matched[i] = matched;
        }
    }
    ctx->merge_global();
    return nullptr;
}

void validate_options(const Options &options)
{
    if(options.source.empty()) throw std::invalid_argument("--source is required");
    if(options.target.empty() && options.query != Query::Within)
        throw std::invalid_argument("--target is required for intersect and intersection");
    if(options.target_object < -1)
        throw std::invalid_argument("target_object must be -1 (all) or a nonnegative source object ID");
    if(options.threads < 0 || options.vertices_per_raster <= 0 ||
       options.bitwidth < 2 || options.bitwidth > 8 || options.distance_km < 0 ||
       options.lowest_layer < 1 || options.unroll_size <= 0 ||
       !std::isfinite(options.confidence) || options.confidence < 0 || options.confidence > 1 ||
       !std::isfinite(options.merge_threshold) || options.merge_threshold < 0 || options.merge_threshold > 1)
        throw std::invalid_argument("Invalid FARM parameters; use tdbase farm --help");
    if(options.backend == Backend::GPU && !gpu_available())
        throw std::runtime_error("This build has no CUDA support; use --cpu or rebuild with -DUSE_GPU=ON");
}

std::size_t convert(const std::string &input, const std::string &output, bool filter)
{
    if(input.empty() || output.empty()) throw std::invalid_argument("Input and output paths are required");
    if(std::filesystem::weakly_canonical(input) == std::filesystem::weakly_canonical(output))
        throw std::invalid_argument("Input and output must be different files");
    std::ifstream check(input);
    if(!check) throw std::runtime_error("Cannot open WKT file: " + input);
    std::string line;
    std::size_t row = 0;
    while(std::getline(check, line)){
        ++row;
        validate_wkt(line, row);
    }
    auto objects = load_polygon_wkt(input.c_str());
    std::vector<Farm *> kept;
    for(Farm *obj : objects){
        bool valid = true;
        for(int i = 0; i < obj->get_num_vertices(); ++i){
            const auto point = obj->get_boundary()->p[i];
            valid &= std::isfinite(point.x) && std::isfinite(point.y);
            if(filter) valid &= point.x > -180 && point.x < 180 && point.y > -90 && point.y < 90;
        }
        if(valid) kept.push_back(obj);
    }
    std::ofstream probe(output, std::ios::binary);
    if(!probe){
        for(Farm *obj : objects) delete obj;
        throw std::runtime_error("Cannot write output: " + output);
    }
    probe.close();
    dump_polygons_to_file(kept, output.c_str());
    for(Farm *obj : objects) delete obj;
    return kept.size();
}
} // namespace

bool gpu_available()
{
#ifdef USE_GPU
    return true;
#else
    return false;
#endif
}

Result execute(const Options &options)
{
    validate_options(options);
    Result result;
    const auto load_start = std::chrono::steady_clock::now();
    const auto sources = validate_binary(options.source);
    if(options.target_object >= 0 && static_cast<std::uint64_t>(options.target_object) >= sources)
        throw std::out_of_range("Source object ID " + std::to_string(options.target_object) +
            " is out of range (source contains " + std::to_string(sources) + " objects)");
    const auto targets = options.target.empty() ? sources : validate_binary(options.target);
    if(!sources || !targets) return {};
    query_context ctx;
    ctx.use_gpu = options.backend == Backend::GPU;
    ctx.collect_results = options.collect_results;
    ctx.source_path = options.source;
    ctx.target_path = options.target;
    // A selected self join loads the source object once and the full dataset as targets.
    if(options.target_object >= 0 && options.target.empty()) ctx.target_path = options.source;
    ctx.query_type = options.query == Query::Within ? QueryType::within :
        (options.query == Query::Intersection ? QueryType::intersection : QueryType::intersect);
    if(options.threads) ctx.num_threads = options.threads;
    ctx.use_hierarchy = options.hierarchy;
    ctx.use_approximation = options.approximation;
    ctx.vpr = options.vertices_per_raster;
    ctx.bitwidth = options.bitwidth;
    ctx.within_distance = options.distance_km;
    ctx.batch_size = options.batch_size;
    ctx.merge_threshold = options.merge_threshold;
    ctx.approx_confidence = options.confidence;
    ctx.NLow = options.lowest_layer;
    ctx.unroll_size = options.unroll_size;
    if(options.target_object >= 0){
        auto source = load_source_object(ctx.source_path, sources, static_cast<std::size_t>(options.target_object));
        ctx.source_objects.push_back(source.get());
        source.release();
    }else{
        ctx.source_objects = load_binary_file(ctx.source_path.c_str(), ctx);
    }
    if(!ctx.target_path.empty()) ctx.target_objects = load_binary_file(ctx.target_path.c_str(), ctx);
    for(const auto *objects : {&ctx.source_objects, &ctx.target_objects}){
        for(Farm *obj : *objects){
            if(!obj->get_holes().empty())
                throw std::invalid_argument("FARM query kernels require polygons without interior rings");
        }
    }
    result.load_milliseconds = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - load_start).count();
    const auto index_start = std::chrono::steady_clock::now();
    CandidateIndex index;
    index.selected_source = options.target_object;
    index.self_join = options.target.empty();
    for(Farm *obj : ctx.source_objects) index.tree.Insert(obj->getMBB()->low, obj->getMBB()->high, obj);
    ctx.target_num = ctx.target_path.empty() ? ctx.source_objects.size() : ctx.target_objects.size();
    ctx.index = 0;
    ctx.run_worker_threads(collect_candidates, nullptr, &index);
    result.index_milliseconds = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - index_start).count();
    result.candidates = ctx.num_pairs;
    if(!ctx.num_pairs) return result;
    std::sort(ctx.object_pairs.begin(), ctx.object_pairs.end());
    if(ctx.collect_results) ctx.matched.resize(ctx.num_pairs);
    if(options.query == Query::Intersection) ctx.areas = new double[ctx.num_pairs]();
    const auto rasterization_start = std::chrono::steady_clock::now();
    preprocess(&ctx);
    result.rasterization_milliseconds = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - rasterization_start).count();
    auto compute_start = std::chrono::steady_clock::now();
    if(ctx.use_gpu){
#ifdef USE_GPU
        preprocess_for_gpu(&ctx);
        if(!ctx.batch_size) ctx.batch_size = std::min(ctx.num_pairs, std::size_t(100000));
        compute_start = std::chrono::steady_clock::now();
        for(ctx.index = 0; ctx.index < ctx.num_pairs; ctx.index += ctx.batch_size){
            ctx.index_end = std::min(ctx.index + ctx.batch_size, ctx.num_pairs);
            if(options.query == Query::Within) cuda_within_polygon(&ctx);
            else if(options.query == Query::Intersection) cuda_intersection(&ctx);
            else cuda_intersect(&ctx);
        }
#endif
    }else{
        ctx.index = 0;
        ctx.target_num = ctx.num_pairs;
        ctx.run_worker_threads(compute_cpu);
    }
    // GPU batch functions copy their results to the host before returning.
    result.compute_milliseconds = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - compute_start).count();
    if(options.query == Query::Intersection){
        ctx.found = 0;
        for(std::size_t i = 0; i < ctx.num_pairs; ++i){
            // FARM's internal shoelace accumulators store twice the area.
            // Normalize once at the public API boundary for every backend.
            ctx.areas[i] *= 0.5;
            ctx.found += ctx.areas[i] > 0;
        }
    }
    result.count = ctx.found;
    if(ctx.collect_results){
        for(std::size_t i = 0; i < ctx.num_pairs; ++i){
            // Area output includes every candidate, including zero-area pairs.
            if(options.query != Query::Intersection && !ctx.matched[i]) continue;
            const auto pair = ctx.object_pairs[i];
            const auto target = pair.second -
                (ctx.use_gpu && !ctx.target_path.empty() ? ctx.source_objects.size() : 0);
            const auto source = options.target_object >= 0
                ? static_cast<uint32_t>(options.target_object) : pair.first;
            result.matches.push_back({source, static_cast<uint32_t>(target), ctx.areas ? ctx.areas[i] : 0});
        }
    }
    ctx.profiler.print();
    return result;
}

std::size_t convert_wkt(const std::string &input, const std::string &output) { return convert(input, output, false); }
std::size_t filter_wkt(const std::string &input, const std::string &output) { return convert(input, output, true); }

} // namespace tdbase::farm2d
