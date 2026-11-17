#include "api.h"
#include <chrono>
#include <boost/program_options.hpp>
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace tdbase::farm2d {
namespace {
void help()
{
    std::cout << "TDBase FARM 2D\n"
        "  tdbase farm convert INPUT.wkt OUTPUT.dt\n"
        "  tdbase farm filter INPUT.wkt OUTPUT.dt\n"
        "  tdbase farm intersect -s SOURCE.dt -t TARGET.dt [options]\n"
        "  tdbase farm intersection -s SOURCE.dt -t TARGET.dt [options]\n"
        "  tdbase farm within -s SOURCE.dt [-t TARGET.dt] -d KM [options]\n"
        "  tdbase join --dim 2 -q intersect|intersection|within --tile1 FILE [--tile2 FILE]\n\n"
        "Options:\n"
        "  -g, --gpu                  Use CUDA (default: CPU)\n"
        "      --cpu                  Use CPU explicitly\n"
        "      --target-object ID     Query one source object (zero-based; default: all)\n"
        "  -h, --hierarchy            Enable hierarchical rasters\n"
        "  -a, --approximation        Approximate query\n"
        "  -n, --threads N            CPU worker threads\n"
        "  -v, --vpr N                Vertices per raster (default: 10)\n"
        "  -w, --bitwidth N           Fill-ratio bits (2..8; default: 3)\n"
        "  -d, --within_distance KM   Nonnegative integer kilometers (default: 10)\n"
        "  -b, --batch_size N         GPU candidate pairs per batch (default: up to 100000)\n"
        "  -m, --merge_threshold F    GPU merge threshold [0,1]\n"
        "  -c, --approx_confidence F  Approximation confidence [0,1]\n"
        "  -l, --NLow N               Lowest hierarchy layer\n"
        "  -u, --unroll_size N        GPU refinement unroll size\n"
        "  -p, --print_result         Print object pairs (plus area for intersection)\n"
        "  -o, --output FILE          Write results to FILE\n"
        "      --format tsv|json      Result format (default: tsv)\n"
        "      --help                Show help (-h means hierarchy)\n\n"
        "IDs are zero-based dataset positions. Within self joins are directed and exclude identity pairs.\n"
        "Intersection exports candidate pairs and areas, including zero areas, in squared input-coordinate units.\n"
        "FARM 2D datasets use a different binary format from 3D tiles.\n"
        "CUDA compiled: " << (gpu_available() ? "yes" : "no") << '\n';
}

void write_result(std::ostream &out, const Result &result, Query query, const std::string &format)
{
    out << std::setprecision(17);
    if(format == "json"){
        out << "{\"count\":" << result.count << ",\"candidates\":" << result.candidates << ",\"matches\":[";
        bool first = true;
        for(const auto &row : result.matches){
            if(!first) out << ',';
            first = false;
            out << "{\"source\":" << row.source << ",\"target\":" << row.target;
            if(query == Query::Intersection) out << ",\"area\":" << row.area;
            out << '}';
        }
        out << "]}\n";
    }else{
        for(const auto &row : result.matches){
            out << row.source << '\t' << row.target;
            if(query == Query::Intersection) out << '\t' << row.area;
            out << '\n';
        }
    }
}
} // namespace

int command(int argc, char **argv)
{
    const auto total_start = std::chrono::steady_clock::now();
    try{
        if(argc < 2){ help(); return 0; }
        std::string action = argv[1];
        if(action == "--help" || action == "help"){ help(); return 0; }
        if(action == "dump_polygons") action = "convert";
        if(action == "filter_polygons") action = "filter";
        if(action == "within_polygon") action = "within";
        if(action == "convert" || action == "filter"){
            if(argc == 3 && std::string(argv[2]) == "--help"){ help(); return 0; }
            if(argc != 4) throw std::invalid_argument("Usage: tdbase farm " + action + " INPUT.wkt OUTPUT.dt");
            const auto count = action == "convert" ? convert_wkt(argv[2], argv[3]) : filter_wkt(argv[2], argv[3]);
            std::cerr << "polygons-written:\t" << count << '\n';
            return 0;
        }
        if(action != "join" && action != "intersect" && action != "intersection" && action != "within")
            throw std::invalid_argument("Unknown FARM command: " + action);
        Options options;
        std::string query_name = action == "join" ? "intersect" : action;
        std::string output, format = "tsv";
        int dimension = 2;
        namespace po = boost::program_options;
        po::options_description desc;
        desc.add_options()
            ("help", "help")("gpu,g", "GPU")("cpu", "CPU")
            ("source,s", po::value<std::string>(&options.source), "source")
            ("target,t", po::value<std::string>(&options.target), "target")
            ("tile1", po::value<std::string>(&options.source), "source")
            ("tile2", po::value<std::string>(&options.target), "target")
            ("target-object", po::value<std::int64_t>(&options.target_object), "source object ID")
            ("query,q", po::value<std::string>(&query_name), "query")
            ("dim", po::value<int>(&dimension), "dimension")
            ("hierarchy,h", "hierarchy")("approximation,a", "approximation")
            ("threads,n", po::value<int>(&options.threads), "threads")
            ("vpr,v", po::value<int>(&options.vertices_per_raster), "vpr")
            ("bitwidth,w", po::value<int>(&options.bitwidth), "bitwidth")
            ("within_distance,d", po::value<int>(&options.distance_km), "distance")
            ("within_dist", po::value<int>(&options.distance_km), "distance")
            ("batch_size,b", po::value<std::size_t>(&options.batch_size), "batch")
            ("merge_threshold,m", po::value<float>(&options.merge_threshold), "merge")
            ("approx_confidence,c", po::value<float>(&options.confidence), "confidence")
            ("NLow,l", po::value<int>(&options.lowest_layer), "lowest layer")
            ("unroll_size,u", po::value<int>(&options.unroll_size), "unroll")
            ("print_result,p", "print results")
            ("output,o", po::value<std::string>(&output), "output")
            ("format", po::value<std::string>(&format), "format");
        po::variables_map vm;
        po::store(po::command_line_parser(argc - 1, argv + 1).options(desc).run(), vm);
        if(vm.count("help")){ help(); return 0; }
        po::notify(vm);
        if(dimension != 2) throw std::invalid_argument("FARM requires --dim 2");
        if(vm.count("gpu") && vm.count("cpu")) throw std::invalid_argument("Choose --gpu or --cpu");
        if(vm.count("threads") && options.threads <= 0) throw std::invalid_argument("--threads must be positive");
        if(vm.count("target-object") && options.target_object < 0)
            throw std::invalid_argument("--target-object must be a nonnegative source object ID");
        if(format != "tsv" && format != "json") throw std::invalid_argument("--format must be tsv or json");
        if(query_name == "within_polygon") query_name = "within";
        if(query_name == "intersect") options.query = Query::Intersect;
        else if(query_name == "intersection") options.query = Query::Intersection;
        else if(query_name == "within") options.query = Query::Within;
        else throw std::invalid_argument("2D query must be intersect, intersection, or within");
        options.backend = vm.count("gpu") ? Backend::GPU : Backend::CPU;
        options.hierarchy = vm.count("hierarchy");
        options.approximation = vm.count("approximation");
        options.collect_results = vm.count("print_result") || !output.empty();
        // Open the result destination before running an expensive query.
        std::ofstream file;
        if(!output.empty()){
            const auto destination = std::filesystem::weakly_canonical(output);
            if(destination == std::filesystem::weakly_canonical(options.source) ||
               (!options.target.empty() && destination == std::filesystem::weakly_canonical(options.target)))
                throw std::invalid_argument("Output must differ from input datasets");
            file.open(output);
            if(!file) throw std::runtime_error("Cannot write output: " + output);
        }
        const auto result = execute(options);
        std::cerr << "candidate-pairs:\t" << result.candidates << "\ncount-found:\t" << result.count << '\n';
        if(vm.count("print_result")) write_result(std::cout, result, options.query, format);
        if(vm.count("print_result")) std::cout.flush();
        if(file.is_open()){
            write_result(file, result, options.query, format);
            file.close();
            if(!file) throw std::runtime_error("Failed writing output: " + output);
        }
        const double total_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - total_start).count();
        std::cerr << std::fixed << std::setprecision(6)
            << "total:\t" << total_seconds << " s\n"
            << "load:\t" << result.load_milliseconds / 1000.0 << " s\n"
            << "index:\t" << result.index_milliseconds / 1000.0 << " s\n"
            << "rasterization:\t" << result.rasterization_milliseconds / 1000.0 << " s\n"
            << "compute:\t" << result.compute_milliseconds / 1000.0 << " s\n";
        return 0;
    }catch(const std::exception &error){
        std::cerr << "FARM: " << error.what() << '\n';
        return 1;
    }
}
} // namespace tdbase::farm2d
