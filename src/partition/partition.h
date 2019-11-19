/*
 * partition.h
 *
 *  Created on: Nov 18, 2019
 *      Author: teng
 */

#ifndef PARTITION_H_
#define PARTITION_H_

#include "../geometry/aab.h"
#include "../PPMC/mymesh.h"
#include "octree.h"
#include "../spatial/spatial.h"

using namespace std;
namespace hispeed{

// partition the space
void partition_space(std::vector<std::string> input_folders,
		std::vector<aab> &tiles, int num_threads, int num_tiles);
void persist_tile(std::vector<aab> &tiles, const char *space_path);
void load_space(std::vector<aab> &tiles, const char *space_path);
void partition_data(std::vector<aab> &tiles, const char *output_folder);


}
#endif /* PARTITION_H_ */
