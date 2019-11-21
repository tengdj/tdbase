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
	std::vector<string> input_folders;
	input_folders.push_back(argv[1]);

	int num_tiles = hispeed::get_num_threads();
	if(argc>2){
		num_tiles = atoi(argv[2]);
	}

	int sample_rate = 100;
	if(argc>3){
		sample_rate = atoi(argv[3]);
	}

	std::vector<aab *> mbbs;
	hispeed::get_mbbs(input_folders, mbbs, hispeed::get_num_threads(), sample_rate);

	SPNode *sp = NULL;
	if(true){
		sp = hispeed::build_sort_partition(mbbs, num_tiles);
		sp->persist("tmp.off");
	}else{
		sp = new SPNode();
		sp->load("tmp.off");
	}
	sp->genTiles(tiles);
	hispeed::persist_tile(tiles, "tiles_");
	tiles.clear();
	for(aab *b:mbbs){
		delete b;
	}
	mbbs.clear();
	delete sp;

	return 0;
}



