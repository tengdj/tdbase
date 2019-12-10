/*
 * joiner.cpp
 *
 *  Created on: Nov 20, 2019
 *      Author: teng
 */

#include <boost/program_options.hpp>
#include <vector>

#include "../join/SpatialJoin.h"
#include "../spatial/himesh.h"
#include "../index/index.h"

using namespace std;
using namespace hispeed;
namespace po = boost::program_options;

int main(int argc, char **argv){
	struct timeval start = get_cur_time();

	string tile1_path("nuclei_tmp.dt");
	string tile2_path("nuclei_tmp.dt");
	bool intersect = false;
	int num_threads = hispeed::get_num_threads();
	size_t max_objects = LONG_MAX;

	po::options_description desc("joiner usage");
	desc.add_options()
		("help,h", "produce help message")
		("gpu,g", "compute with GPU")
		("intersect,i", "do intersection instead of join")
		("tile1", po::value<string>(&tile1_path), "path to tile 1")
		("tile2", po::value<string>(&tile2_path), "path to tile 2")
		("threads,n", po::value<int>(&num_threads), "number of threads")
		("max_objects,m", po::value<size_t>(&max_objects), "max number of objects in a tile")
		;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	if (vm.count("help")) {
		cout << desc << "\n";
		return 0;
	}
	po::notify(vm);

	geometry_computer *gc = new geometry_computer();
	if(vm.count("gpu")){
		initialize();
		gc->init_gpus();
	}
	if(vm.count("intersect")){
		intersect = true;
	}
	if(vm.count("threads")&&num_threads>0){
		gc->set_thread_num(num_threads);
	}

	Tile *tile1 = new Tile(tile1_path.c_str(), max_objects);
	Tile *tile2 = tile1;
	if(vm.count("tile2")&&tile1_path!=tile2_path){
		tile2 = new Tile(tile2_path.c_str(), max_objects);
	}
	assert(tile1&&tile2);
	logt("load tiles", start);

	SpatialJoin *joiner = new SpatialJoin(gc);
	if(intersect){
		joiner->intersect(tile1, tile2);
	}else{
		vector<pair<Tile *, Tile *>> tile_pairs;
		for(int i=0;i<8;i++){
			tile_pairs.push_back(pair<Tile *, Tile *>(tile1, tile2));
		}
		joiner->nearest_neighbor_batch(tile_pairs, 8);
	}
	logt("join", start);

	delete tile1;
	if(tile2!=tile1){
		delete tile2;
	}
	delete joiner;
	delete gc;
	logt("cleaning", start);

}
