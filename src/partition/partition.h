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
void get_voxels(std::vector<std::string> &input_folders,
		std::vector<weighted_aab *> &voxels, const int num_threads, const int sample_rate);


}
#endif /* PARTITION_H_ */
