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

	query_context ctx;
	int base_lod = 0;
	int lod_gap = 50;
	int top_lod = 100;

	po::options_description desc("joiner usage");
	desc.add_options()
		("help,h", "produce help message")
		("query,q", po::value<string>(&ctx.query_type),"query type can be intersect|nn|within")
		("knn", po::value<int>(&ctx.knn), "the K value for NN query")
		("tile1", po::value<string>(&tile1_path), "path to tile 1")
		("tile2", po::value<string>(&tile2_path), "path to tile 2")
		("cn", po::value<int>(&ctx.num_compute_thread), "number of threads for geometric computation for each tile")
		("threads,n", po::value<int>(&ctx.num_thread), "number of threads for processing tiles")
		("repeat,r", po::value<int>(&ctx.repeated_times), "repeat tiles")
		("max_objects1", po::value<size_t>(&ctx.max_num_objects1), "max number of objects in tile 1")
		("max_objects2", po::value<size_t>(&ctx.max_num_objects2), "max number of objects in tile 2")
		("base_lod", po::value<int>(&base_lod), "the base lod for progressive decoding polyhedral")
		("top_lod", po::value<int>(&top_lod), "the top lod for progressive decoding polyhedral")
		("lod_gap", po::value<int>(&lod_gap), "the lod gap for progressive decoding polyhedral")
		("lod", po::value<std::vector<std::string>>()->multitoken()->
		        zero_tokens()->composing(), "the lods need be processed")
		("aabb", "calculate distance with aabb")
		("gpu,g", "compute with GPU")
		("multiple_mbb,m", "using shape-aware indexing with multiple MBB")
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
	if(vm.count("multiple_mbb")){
		ctx.use_multimbb = true;
	}
	if(ctx.num_compute_thread>0){
		gc->set_thread_num(ctx.num_compute_thread);
	}

	if(ctx.query_type!="intersect"&&ctx.query_type!="nn"&&ctx.query_type!="within"){
		cout <<"error query type: "<< ctx.query_type <<endl;
		return 0;
	}


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
		Tile *tile1, *tile2;
		if(vm.count("tile2")){
			tile1 = new Tile(tile1_path.c_str(), ctx.max_num_objects1);
			tile2 = new Tile(tile2_path.c_str(), ctx.max_num_objects2);
		}else{
			tile1 = new Tile(tile1_path.c_str(), LONG_MAX);
			tile2 = tile1;
		}

		assert(tile1&&tile2);
		if(!ctx.use_multimbb){
			tile1->disable_innerpart();
			tile2->disable_innerpart();
		}
		tile_pairs.push_back(pair<Tile *, Tile *>(tile1, tile2));
	}
	logt("load tiles", start);

	SpatialJoin *joiner = new SpatialJoin(gc,ctx);
	joiner->join(tile_pairs, ctx);

	double join_time = hispeed::get_time_elapsed(start,false);
	logt("join", start);
	tile_pairs.clear();
	joiner->report_time(join_time);
	delete joiner;
	delete gc;
	logt("cleaning", start);

}
