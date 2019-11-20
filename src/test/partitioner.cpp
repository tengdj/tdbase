/*
 * partitioner.cpp
 *
 *  Created on: Nov 18, 2019
 *      Author: teng
 */

#include "../partition/partition.h"
#include <boost/program_options.hpp>

using namespace hispeed;
using namespace std;

// main method
int main(int argc, char** argv) {
	std::vector<aab> tiles;

//	// parse mbb
//	std::vector<string> input_folders;
//	input_folders.push_back(argv[1]);
//	std::vector<aab *> mbbs;
//	hispeed::get_mbbs(input_folders, mbbs, hispeed::get_num_threads());
//
//	int num_tiles = hispeed::get_num_threads();
//	if(argc>2){
//		num_tiles = atoi(argv[2]);
//	}
//
//	SPNode *sp = hispeed::build_sort_partition(mbbs, num_tiles);
//	sp->persist("tmp");
//	sp->genTiles(tiles);
//	delete sp;
	SPNode *sp = new SPNode();
	sp->load("tmp");
	sp->genTiles(tiles);
	hispeed::persist_tile(tiles, "tiles_");
	tiles.clear();
//	for(aab *b:mbbs){
//		delete b;
//	}
//	mbbs.clear();

	return 0;
}



