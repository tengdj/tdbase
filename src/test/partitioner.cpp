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
	if(argc<2){
		cerr<<"usage: parition path/to/folder [num_threads] [sample_rate]"<<endl;
		exit(0);
	}
	std::vector<aab> tiles;

	// parse mbb
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

	struct timeval start = get_cur_time();

	SPNode *sp = NULL;
	const char *partition_path = "tmp.part";
	if(!hispeed::file_exist(partition_path)){
		std::vector<Voxel *> voxels;
		hispeed::get_voxels(input_folders, voxels,
				hispeed::get_num_threads(), sample_rate);
		cout<<"getting voxels of objects takes "<<get_time_elapsed(start, true)<<" ms"<<endl;
		sp = hispeed::build_sort_partition(voxels, num_tiles);
		cout<<"generating partitions takes "<<get_time_elapsed(start, true)<<" ms"<<endl;
		for(Voxel *v:voxels){
			delete v;
		}
		voxels.clear();
		sp->persist(partition_path);
		cout<<"persist partitions takes "<<get_time_elapsed(start, true)<<" ms"<<endl;
	}else{
		sp = new SPNode();
		sp->load(partition_path);
		cout<<"loading partitions takes "<<get_time_elapsed(start, true)<<" ms"<<endl;
	}
	sp->genTiles(tiles);
	hispeed::persist_tile(tiles, "tiles_");
	tiles.clear();

	delete sp;
	return 0;
}



