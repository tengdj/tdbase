/*
 * partition.h
 *
 *  Created on: Nov 18, 2019
 *      Author: teng
 */

#ifndef PARTITION_H_
#define PARTITION_H_
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <map>
#include <vector>
#include <cstdlib>
#include <algorithm>

#include "../geometry/aab.h"
#include "../PPMC/mymesh.h"
#include "../spatial/spatial.h"
#include "../index/index.h"

using namespace std;
namespace hispeed{

// partition the space
void get_mbbs(std::vector<std::string> &input_folders,
		std::vector<aab *> &mbbs, const int num_threads, const int sample_rate);

void persist_tile(std::vector<aab> &tiles, const char *space_path);
void load_space(std::vector<aab> &tiles, const char *space_path);
void partition_data(std::vector<aab> &tiles, const char *output_folder);

}
#endif /* PARTITION_H_ */
