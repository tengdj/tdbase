/*
 * partitioner.cpp
 *
 *  Created on: Nov 18, 2019
 *      Author: teng
 */

#include "../partition/partition.h"

using namespace hispeed;
using namespace std;

// main method
int main(int argc, char** argv) {

	std::vector<string> input_folders;
	//input_folders.push_back("/home/teng/project/HiSPEED/build/nuclei");
	//input_folders.push_back("/home/teng/project/HiSPEED/build/vessels");
	input_folders.push_back(argv[1]);
	std::vector<aab *> mbbs;
	hispeed::get_mbbs(input_folders, mbbs, hispeed::get_num_threads());
	std::vector<aab> tiles;

//	OctreeNode *octree = hispeed::build_ctree(mbbs, 100);
//	octree->genTiles(tiles);
//	delete octree;
	int num_tiles = 100;
	if(argc>2){
		num_tiles = atoi(argv[2]);
	}

	SPNode *sp = hispeed::build_sort_partition(mbbs, num_tiles);
	sp->genTiles(tiles);

	hispeed::persist_tile(tiles, "tiles");
	tiles.clear();
	for(aab *b:mbbs){
		delete b;
	}
	mbbs.clear();
	return 0;
}



