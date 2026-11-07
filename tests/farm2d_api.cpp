#include "farm2d/api.h"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

using namespace tdbase::farm2d;

static void require(bool condition, const char *message)
{
    if(!condition) throw std::runtime_error(message);
}

int main(int argc, char **argv)
{
    try{
        const std::filesystem::path work = argv[1];
        std::filesystem::create_directories(work);
        std::ofstream(work / "a.wkt") << "POLYGON ((0 0,1 0,1 1,0 1,0 0))\n";
        std::ofstream(work / "b.wkt") << "POLYGON ((0.5 0.5,1.5 0.5,1.5 1.5,0.5 1.5,0.5 0.5))\n";
        std::ofstream(work / "near.wkt") << "POLYGON ((1.01 0,2.01 0,2.01 1,1.01 1,1.01 0))\n";
        require(convert_wkt((work / "a.wkt").string(), (work / "a.dt").string()) == 1, "convert a");
        require(convert_wkt((work / "b.wkt").string(), (work / "b.dt").string()) == 1, "repeat conversion");
        require(convert_wkt((work / "near.wkt").string(), (work / "near.dt").string()) == 1, "convert near");
        Options options;
        options.source = (work / "a.dt").string();
        options.target = (work / "b.dt").string();
        options.threads = 2;
        options.vertices_per_raster = 1;
        for(Backend backend : {Backend::CPU, Backend::GPU}){
            if(backend == Backend::GPU && !gpu_available()) continue;
            options.backend = backend;
            for(bool hierarchy : {false, true}){
                options.hierarchy = hierarchy;
                for(int bits : {2, 4, 8}){
                    options.bitwidth = bits;
                    options.query = Query::Intersect;
                    auto intersection = execute(options);
                    require(intersection.count == 1 && intersection.matches.size() == 1, "intersect API");
                    require(intersection.matches[0].source == 0 && intersection.matches[0].target == 0, "local IDs");
                    options.query = Query::Intersection;
                    auto area = execute(options);
                    require(area.matches.size() == 1 && std::abs(area.matches[0].area - .25) < 1e-6, "numeric intersection area");
                    options.query = Query::Within;
                    options.target = (work / "near.dt").string();
                    options.distance_km = 1;
                    require(execute(options).count == 0, "within rejects distant pair");
                    options.distance_km = 2;
                    auto within = execute(options);
                    require(within.count == 1 && within.matches.size() == 1, "within API output");
                    options.target = (work / "b.dt").string();
                }
            }
        }
        std::cout << "API repeat calls, CPU/GPU, flat/hierarchy, bitwidths, areas and IDs passed\n";
        return 0;
    }catch(const std::exception &error){
        std::cerr << error.what() << '\n';
        return 1;
    }
}
