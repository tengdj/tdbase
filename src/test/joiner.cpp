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

	query_context ctx;
	size_t max_objects = LONG_MAX;
	int base_lod = 0;
	int lod_gap = 50;
	int top_lod = 100;


	po::options_description desc("joiner usage");
	desc.add_options()
		("help,h", "produce help message")
		("gpu,g", "compute with GPU")
		("query,q", po::value<string>(&query),"query type can be intersect|nn|within")
		("tile1", po::value<string>(&tile1_path), "path to tile 1")
		("tile2", po::value<string>(&tile2_path), "path to tile 2")
		("threads,n", po::value<int>(&ctx.num_thread), "number of threads")
		("rn", po::value<int>(&ctx.num_repeated_thread), "number of threads for repeating jobs")
		("repeat,r", po::value<int>(&ctx.repeated_times), "repeat tiles")
		("max_objects,m", po::value<size_t>(&max_objects), "max number of objects in a tile")
		("base_lod", po::value<int>(&base_lod), "the base lod for progressive decoding polyhedral")
		("top_lod", po::value<int>(&top_lod), "the top lod for progressive decoding polyhedral")
		("lod_gap", po::value<int>(&lod_gap), "the lod gap for progressive decoding polyhedral")
		("lod", po::value<std::vector<std::string>>()->multitoken()->
		        zero_tokens()->composing(), "the lods need be processed")
		("aabb", "calculate distance with aabb")
		("max_dist", po::value<double>(&ctx.max_dist), "the maximum distance for within query")
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
	if(vm.count("aabb")){
		ctx.use_aabb = true;
	}
	if(vm.count("threads")&&ctx.num_thread>0){
		gc->set_thread_num(ctx.num_thread);
	}


	SpatialJoin *joiner = new SpatialJoin(gc);
	if(vm.count("lod")){
		for(string l:vm["lod"].as<std::vector<std::string>>()){
			ctx.lods.push_back(atoi(l.c_str()));
		}
	}else{
		for(int l=base_lod;l<=top_lod;l+=lod_gap){
			ctx.lods.push_back(l);
		}
		if(ctx.lods[ctx.lods.size()-1]<top_lod){
			ctx.lods.push_back(top_lod);
		}
	}

	vector<pair<Tile *, Tile *>> tile_pairs;
	for(int i=0;i<ctx.repeated_times;i++){
		Tile *tile1 = new Tile(tile1_path.c_str(), max_objects);
		Tile *tile2 = tile1;
		if(vm.count("tile2")){
			tile2 = new Tile(tile2_path.c_str(), max_objects);
		}
		assert(tile1&&tile2);
		tile_pairs.push_back(pair<Tile *, Tile *>(tile1, tile2));
		tile2->retrieve_all();
		tile2->advance_all(100);
		for(int i=0;i<tile2->num_objects();i++){
			tile2->get_mesh_wrapper(i)->writeMeshOff();
		}

	}
	logt("load tiles", start);

	if(query=="intersect"){
		joiner->intersect_batch(tile_pairs, ctx);
	}else if(query=="nn"){
		joiner->nearest_neighbor_batch(tile_pairs, ctx);
	}else if(query=="within"){
		joiner->within_batch(tile_pairs, ctx);
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
