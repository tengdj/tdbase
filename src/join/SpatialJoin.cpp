/*
 * SpatialJoin.cpp
 *
 *  Created on: Nov 11, 2019
 *      Author: teng
 */

#include <math.h>
#include <map>
#include <tuple>
#include <string.h>
#include "SpatialJoin.h"

using namespace std;

namespace tdbase{


/*facilitate functions*/

size_t get_pair_num(vector<candidate_entry *> &candidates){
	size_t pair_num = 0;
	for(candidate_entry *p:candidates){
		for(candidate_info &c:p->candidates){
			pair_num += c.voxel_pairs.size();
		}
	}
	return pair_num;
}

size_t get_candidate_num(vector<candidate_entry *> &candidates){
	size_t candidate_num = 0;
	for(candidate_entry *p:candidates){
		candidate_num += p->candidates.size();
	}
	return candidate_num;
}

SpatialJoin::SpatialJoin(){
	computer = new geometry_computer();
	if(config.use_gpu){
#ifdef USE_GPU
		initialize();
		computer->init_gpus();
#endif
	}
	if(config.num_compute_thread>0){
		computer->set_thread_num(config.num_compute_thread);
	}
}

SpatialJoin::~SpatialJoin(){
	delete computer;
}

void SpatialJoin::decode_data(query_context &ctx){
	struct timeval start = tdbase::get_cur_time();
	// decode the objects to current lod
	for(candidate_entry *c:ctx.candidates){
		HiMesh_Wrapper *wrapper1 = c->mesh_wrapper;
		wrapper1->decode_to(ctx.cur_lod);
		for(candidate_info &info:c->candidates){
			info.mesh_wrapper->decode_to(ctx.cur_lod);
		}// end for distance_candiate list
	}// end for candidates
	ctx.decode_time += logt("decode data", start);
}

void SpatialJoin::packing_data(query_context &ctx){
	struct timeval start = tdbase::get_cur_time();

	// initialize the result space
	ctx.gp.pair_num = get_pair_num(ctx.candidates);
	if(ctx.gp.results){
		delete []ctx.gp.results;
	}
	ctx.gp.results = new result_container[ctx.gp.pair_num];
	for(int i=0;i<ctx.gp.pair_num;i++){
		ctx.gp.results[i].distance = 0;
		ctx.gp.results[i].intersected = false;
	}

	if(config.use_aabb){
		// build the AABB tree
		for(candidate_entry *c:ctx.candidates){
			for(candidate_info &info:c->candidates){
				info.mesh_wrapper->get_mesh()->get_aabb_tree_triangle();
			}
			c->mesh_wrapper->get_mesh()->get_aabb_tree_triangle();
		}
		ctx.packing_time += logt("building aabb tree", start);
	}else{
		// for progressive query, get the triangles
		ctx.gp.element_num = 0;
		ctx.gp.element_pair_num = 0;
		map<Voxel *, uint32_t> voxel_offset_map;

		for(candidate_entry *c:ctx.candidates){
			HiMesh_Wrapper *wrapper1 = c->mesh_wrapper;
			for(candidate_info &info:c->candidates){
				for(voxel_pair &vp:info.voxel_pairs){
					ctx.gp.element_pair_num += vp.v1->num_triangles*vp.v2->num_triangles;
					// update the voxel offset map
					for(int i=0;i<2;i++){
						Voxel *tv = (i==0?vp.v1:vp.v2);
						// first time visited this voxel
						if(voxel_offset_map.find(tv)==voxel_offset_map.end()){
							// the voxel is inserted into the map with an offset
							voxel_offset_map[tv] = ctx.gp.element_num;
							ctx.gp.element_num += tv->num_triangles;
						}
					}
				}
			}
		}

		ctx.gp.allocate_buffer();

		// now we allocate the space and store the data in the buffer
		for (map<Voxel *, uint32_t>::iterator it=voxel_offset_map.begin(); it!=voxel_offset_map.end(); ++it){
			Voxel *v = it->first;
			if(v->num_triangles > 0){
				memcpy(ctx.gp.data+voxel_offset_map[v]*9, v->triangles, v->num_triangles*9*sizeof(float));
				memcpy(ctx.gp.hausdorff+voxel_offset_map[v]*2, v->hausdorff, v->num_triangles*2*sizeof(float));
			}
		}

		// organize the data for computing
		int index = 0;
		for(candidate_entry *c:ctx.candidates){
			for(candidate_info &info:c->candidates){
				for(voxel_pair &vp:info.voxel_pairs){
					ctx.gp.offset_size[4*index] = voxel_offset_map[vp.v1];
					ctx.gp.offset_size[4*index+1] = vp.v1->num_triangles;
					ctx.gp.offset_size[4*index+2] = voxel_offset_map[vp.v2];
					ctx.gp.offset_size[4*index+3] = vp.v2->num_triangles;
					index++;
				}
			}
		}
		assert(index==ctx.gp.pair_num);
		voxel_offset_map.clear();
		ctx.packing_time += logt("organizing data with %ld elements and %ld element pairs", start, ctx.gp.element_num, ctx.gp.element_pair_num);
	}
}


/*
 * the main function for conducting the join operation
 *
 * */
void SpatialJoin::join(Tile *tile1, Tile *tile2){
	struct timeval very_start = get_cur_time();
	query_context ctx;
	ctx.obj_count = tile1->num_objects();

	// filtering with MBBs to get the candidate list
	index_retrieval(tile1, tile2, ctx);

	// now we start to conduct query with progressive level of details
	for(uint32_t lod:config.lods){
		ctx.cur_lod = lod;
		struct timeval iter_start = get_cur_time();

		// print some statistics for this round of evaluations
		size_t candidate_num = get_candidate_num(ctx.candidates);
		const int pair_num = get_pair_num(ctx.candidates);
		if(pair_num==0){
			break;
		}
		log("%ld polyhedron has %d candidates %d voxel pairs %.2f voxel pairs per candidate",
				ctx.candidates.size(), candidate_num, pair_num, (1.0*pair_num)/ctx.candidates.size());

		// decode the corresponding polyhedrons into current LOD
		decode_data(ctx);

		// prepare data for conducting geometric computation (AABBtree building included)
		packing_data(ctx);

		// truly conduct the geometric computations
		geometric_computation(ctx);

		// now update the candidate list with the latest information calculated with this LOD
		evaluate_candidate_lists(ctx);

		logt("evaluating with lod %d", iter_start, lod);
		log("");
		if(lod==config.highest_lod()){
			assert(ctx.candidates.size()==0);
		}
	}

	ctx.overall_time = get_time_elapsed(very_start, false);
	if(config.print_result){
		ctx.print_result();
	}
	ctx.report();
}

}


