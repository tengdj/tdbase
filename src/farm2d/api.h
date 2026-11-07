#ifndef TDBASE_FARM2D_API_H
#define TDBASE_FARM2D_API_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace tdbase::farm2d {

enum class Query { Intersect, Intersection, Within };
enum class Backend { CPU, GPU };

struct Options {
    Query query = Query::Intersect;
    Backend backend = Backend::CPU;
    std::string source;
    std::string target; // Empty means a directed self join for Within.
    int threads = 0; // 0 selects the available CPU count.
    bool hierarchy = false;
    bool approximation = false;
    int vertices_per_raster = 10;
    int bitwidth = 4;
    int distance_km = 10;
    std::size_t batch_size = 0;
    float merge_threshold = 0.9f;
    float confidence = 0.9f;
    int lowest_layer = 1;
    int unroll_size = 16;
    bool collect_results = true;
};

struct Match {
    std::uint32_t source;
    std::uint32_t target;
    double area = 0; // Intersection only; units are squared input coordinates.
};

struct Result {
    std::size_t candidates = 0;
    std::size_t count = 0;
    std::vector<Match> matches;
};

bool gpu_available(); // Whether this build includes FARM CUDA support.
Result execute(const Options &options);
std::size_t convert_wkt(const std::string &input, const std::string &output);
std::size_t filter_wkt(const std::string &input, const std::string &output);
int command(int argc, char **argv);

} // namespace tdbase::farm2d
#endif
