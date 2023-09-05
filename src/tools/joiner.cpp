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

	char path1[256];
	char path2[256];

	vector<pair<Tile *, Tile *>> tile_pairs;
	for(int i=0;i<global_ctx.repeated_times;i++){
		Tile *tile1, *tile2;
		if(global_ctx.tile2_path.size()>0){
			if(global_ctx.use_raw){
				sprintf(path1, "%s.raw", global_ctx.tile1_path.c_str());
				sprintf(path2, "%s.raw", global_ctx.tile2_path.c_str());
			}else {
				sprintf(path1, "%s", global_ctx.tile1_path.c_str());
				sprintf(path2, "%s", global_ctx.tile2_path.c_str());
			}
			tile1 = new Tile(path1, global_ctx.max_num_objects1, global_ctx.use_raw?RAW:COMPRESSED, false);
			tile2 = new Tile(path2, global_ctx.max_num_objects2, global_ctx.use_raw?RAW:COMPRESSED, false);
		}else{
			if(global_ctx.use_raw){
				sprintf(path1, "%s.raw", global_ctx.tile1_path.c_str());
			}else {
				sprintf(path1, "%s", global_ctx.tile1_path.c_str());
			}
			tile1 = new Tile(path1, LONG_MAX, global_ctx.use_raw?RAW:COMPRESSED, false);
			tile2 = tile1;
		}
		assert(tile1&&tile2);
		tile_pairs.push_back(pair<Tile *, Tile *>(tile1, tile2));
	}
	logt("create tiles", start);

#pragma omp parallel for
	for(int i=0;i<global_ctx.repeated_times;i++){
		Tile *t1 = tile_pairs[i].first;
		Tile *t2 = tile_pairs[i].second;
		t1->init();
		if(t2 != t1){
			t2->init();
		}
	}
	logt("init tiles", start);

	SpatialJoin *joiner = new SpatialJoin(gc);
	joiner->join(tile_pairs);
	double join_time = hispeed::get_time_elapsed(start,false);
	logt("join", start);
	tile_pairs.clear();
	delete joiner;
	delete gc;
}
