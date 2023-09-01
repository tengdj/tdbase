/*
 * joiner.cpp
 *
 *  Created on: Nov 20, 2019
 *      Author: teng
 */

#include <vector>

#include "SpatialJoin.h"
#include "himesh.h"
#include "index.h"

using namespace std;
using namespace hispeed;

int main(int argc, char **argv){
	struct timeval start = get_cur_time();

	global_ctx = parse_args(argc, argv);

	geometry_computer *gc = new geometry_computer();
	if(global_ctx.use_gpu){
#ifdef USE_GPU
		initialize();
		gc->init_gpus();
#endif
	}
	if(global_ctx.num_compute_thread>0){
		gc->set_thread_num(global_ctx.num_compute_thread);
	}

	vector<pair<Tile *, Tile *>> tile_pairs;
	for(int i=0;i<global_ctx.repeated_times;i++){
		Tile *tile1, *tile2;
		if(global_ctx.tile2_path.size()>0){
			tile1 = new Tile(global_ctx.tile1_path.c_str(), global_ctx.max_num_objects1, global_ctx.use_raw?RAW:COMPRESSED);
			tile2 = new Tile(global_ctx.tile2_path.c_str(), global_ctx.max_num_objects2, global_ctx.use_raw?RAW:COMPRESSED);
		}else{
			tile1 = new Tile(global_ctx.tile1_path.c_str(), LONG_MAX, global_ctx.use_raw?RAW:COMPRESSED);
			tile2 = tile1;
		}
		assert(tile1&&tile2);
		tile_pairs.push_back(pair<Tile *, Tile *>(tile1, tile2));
	}
	logt("load tiles", start);

	SpatialJoin *joiner = new SpatialJoin(gc);
	joiner->join(tile_pairs);
	double join_time = hispeed::get_time_elapsed(start,false);
	logt("join", start);
	tile_pairs.clear();
	delete joiner;
	delete gc;
}
