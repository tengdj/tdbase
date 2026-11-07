#pragma once
#include "../api.h"
#include <string>
#include <vector>

inline int farm2d_legacy_main(int argc, char **argv, const char *action, bool gpu = false)
{
    std::vector<std::string> args{argv[0], action};
    if(std::string(action) != "convert" && std::string(action) != "filter")
        args.emplace_back(gpu ? "--gpu" : "--cpu");
    for(int i = 1; i < argc; ++i) args.emplace_back(argv[i]);
    std::vector<char *> pointers;
    for(auto &arg : args) pointers.push_back(arg.data());
    return tdbase::farm2d::command(static_cast<int>(pointers.size()), pointers.data());
}
