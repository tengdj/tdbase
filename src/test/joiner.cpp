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
	string query("intersect");
	bool ispeed = false;
	int num_threads = hispeed::get_num_threads();
	int num_repeat_threads = hispeed::get_num_threads();

	size_t max_objects = LONG_MAX;
	int base_lod = 0;
	int lod_gap = 50;
	int top_lod = 100;
	int repeated = 1;
	double max_dist = 100;

	po::options_description desc("joiner usage");
	desc.add_options()
		("help,h", "produce help message")
		("gpu,g", "compute with GPU")
		("query,q", po::value<string>(&query),"query type can be intersect|nn|within")
		("tile1", po::value<string>(&tile1_path), "path to tile 1")
		("tile2", po::value<string>(&tile2_path), "path to tile 2")
		("threads,n", po::value<int>(&num_threads), "number of threads")
		("rn", po::value<int>(&num_repeat_threads), "number of threads for repeating jobs")
		("repeat,r", po::value<int>(&repeated), "repeat tiles")
		("max_objects,m", po::value<size_t>(&max_objects), "max number of objects in a tile")
		("base_lod", po::value<int>(&base_lod), "the base lod for progressive decoding polyhedral")
		("top_lod", po::value<int>(&top_lod), "the top lod for progressive decoding polyhedral")
		("lod_gap", po::value<int>(&lod_gap), "the lod gap for progressive decoding polyhedral")
		("lod", po::value<std::vector<std::string>>()->multitoken()->
		        zero_tokens()->composing(), "the lods need be processed")
		("ispeed", "run in ispeed mode")
		("max_dist", po::value<double>(&max_dist), "the maximum distance for within query")
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
#ifdef USE_GPU
		initialize();
		gc->init_gpus();
#endif
	}
	if(vm.count("ispeed")){
		ispeed = true;
	}
	if(vm.count("threads")&&num_threads>0){
		gc->set_thread_num(num_threads);
	}


	SpatialJoin *joiner = new SpatialJoin(gc);
	if(vm.count("lod_gap")){
		joiner->set_lod_gap(lod_gap);
	}
	if(vm.count("base_lod")){
		joiner->set_base_lod(base_lod);
	}
	if(vm.count("top_lod")){
		joiner->set_top_lod(top_lod);
	}
	if(vm.count("lod")){
		vector<int> lods;
		for(string l:vm["lod"].as<std::vector<std::string>>()){
			lods.push_back(atoi(l.c_str()));
		}
		joiner->set_lods(lods);
	}

	vector<pair<Tile *, Tile *>> tile_pairs;
	for(int i=0;i<repeated;i++){
		Tile *tile1 = new Tile(tile1_path.c_str(), max_objects);
		Tile *tile2 = tile1;
		if(vm.count("tile2")){
			tile2 = new Tile(tile2_path.c_str(), max_objects);
		}
		assert(tile1&&tile2);
		tile_pairs.push_back(pair<Tile *, Tile *>(tile1, tile2));
	}
	logt("load tiles", start);

	if(query=="intersect"){
		joiner->intersect_batch(tile_pairs, num_repeat_threads);
	}else if(query=="nn"){
		joiner->nearest_neighbor_batch(tile_pairs, num_repeat_threads, ispeed);
	}else if(query=="within"){
		joiner->within_batch(tile_pairs, num_repeat_threads, max_dist);
	}else{
		cout <<"error query type"<<endl<< desc << "\n";
		return 0;
	}
	double join_time = hispeed::get_time_elapsed(start,false);
	logt("join", start);
	tile_pairs.clear();
	joiner->report_time(join_time);
	delete joiner;
	delete gc;
	logt("cleaning", start);

}
