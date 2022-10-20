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

namespace hispeed{


/*facilitate functions*/

size_t get_pair_num(vector<candidate_entry> &candidates){
	size_t pair_num = 0;
	for(candidate_entry &p:candidates){
		for(candidate_info &c:p.candidates){
			pair_num += c.voxel_pairs.size();
		}
	}
	return pair_num;
}

size_t get_candidate_num(vector<candidate_entry> &candidates){
	size_t candidate_num = 0;
	for(candidate_entry &p:candidates){
		candidate_num += p.candidates.size();
	}
	return candidate_num;
}

SpatialJoin::SpatialJoin(geometry_computer *c){
	assert(c);
	computer = c;
}

SpatialJoin::~SpatialJoin(){

}

range SpatialJoin::update_voxel_pair_list(vector<voxel_pair> &voxel_pairs, double minmaxdist){
	range ret;
	ret.mindist = DBL_MAX;
	ret.maxdist = minmaxdist;
	// some voxel pair is farther than this one
	for(auto vp_iter = voxel_pairs.begin();vp_iter!=voxel_pairs.end();){
		// a closer voxel pair already exist
		if(vp_iter->dist.mindist > minmaxdist){
			// evict this unqualified voxel pairs
			voxel_pairs.erase(vp_iter);
		}else{
			ret.mindist = min(ret.mindist, vp_iter->dist.mindist);
			vp_iter++;
		}
	}
	return ret;
}

void SpatialJoin::decode_data(vector<candidate_entry> &candidates, query_context &ctx){
	for(candidate_entry &c:candidates){
		//print_candidate(c);
		HiMesh_Wrapper *wrapper1 = c.mesh_wrapper;
		for(candidate_info &info:c.candidates){
			for(voxel_pair &vp:info.voxel_pairs){
				assert(vp.v1&&vp.v2);
				// not filled yet
				if(vp.v1->data.find(ctx.cur_lod)==vp.v1->data.end()){
					// ensure the mesh is extracted
					ctx.tile1->decode_to(wrapper1->id, ctx.cur_lod);
				}

				if(vp.v2->data.find(ctx.cur_lod)==vp.v2->data.end()){
					ctx.tile2->decode_to(info.mesh_wrapper->id, ctx.cur_lod);
				}
			}// end for voxel_pairs
		}// end for distance_candiate list
	}// end for candidates
}

void SpatialJoin::fill_voxels(vector<candidate_entry> &candidates, query_context &ctx, element_type etype){
	for(candidate_entry &c:candidates){
		HiMesh_Wrapper *wrapper1 = c.mesh_wrapper;
		for(candidate_info &info:c.candidates){
			for(voxel_pair &vp:info.voxel_pairs){
				assert(vp.v1&&vp.v2);
				// not filled yet
				if(vp.v1->data.find(ctx.cur_lod)==vp.v1->data.end()){
					// ensure the mesh is extracted
					wrapper1->fill_voxels(etype);
				}
				if(vp.v2->data.find(ctx.cur_lod)==vp.v2->data.end()){
					info.mesh_wrapper->fill_voxels(etype);
				}
			}// end for voxel_pairs
		}// end for distance_candiate list
	}// end for candidates
}

void SpatialJoin::check_intersection(vector<candidate_entry> &candidates, query_context &ctx){
	struct timeval start = hispeed::get_cur_time();

	const int pair_num = get_pair_num(candidates);

	ctx.results = new result_container[pair_num];
	for(int i=0;i<pair_num;i++){
		ctx.results[i].result.intersected = false;
	}
	decode_data(candidates, ctx);
	ctx.decode_time += logt("decode data", start);

	if(ctx.use_aabb){
		// build the AABB tree
		for(candidate_entry &c:candidates){
			for(candidate_info &info:c.candidates){
				info.mesh_wrapper->mesh->get_aabb_tree_triangle();
			}
		}
		ctx.packing_time += logt("building aabb tree", start);

		int index = 0;
		for(candidate_entry &c:candidates){
			c.mesh_wrapper->mesh->get_segments();
			for(candidate_info &info:c.candidates){
				assert(info.voxel_pairs.size()==1);
				ctx.results[index++].result.intersected = c.mesh_wrapper->mesh->intersect_tree(info.mesh_wrapper->mesh);
			}// end for candidate list
		}// end for candidates
		// clear the trees for current LOD
		for(candidate_entry &c:candidates){
			for(candidate_info &info:c.candidates){
				info.mesh_wrapper->mesh->clear_aabb_tree();
			}
		}
		ctx.computation_time += logt("computation for distance computation", start);
	}else{
		fill_voxels(candidates, ctx, DT_Triangle);
		map<Voxel *, std::pair<uint, uint>> voxel_map;
		size_t element_num = 0;
		size_t element_pair_num = 0;

		for(candidate_entry &c:candidates){
			HiMesh_Wrapper *wrapper1 = c.mesh_wrapper;
			for(candidate_info &info:c.candidates){
				for(voxel_pair &vp:info.voxel_pairs){
					assert(vp.v1&&vp.v2);
					element_pair_num += vp.v1->size[ctx.cur_lod]*vp.v2->size[ctx.cur_lod];
					// update the voxel map
					for(int i=0;i<2;i++){
						Voxel *tv = i==0?vp.v1:vp.v2;
						if(voxel_map.find(tv)==voxel_map.end()){
							std::pair<uint, uint> p;
							p.first = element_num;
							p.second = tv->size[ctx.cur_lod];
							element_num += tv->size[ctx.cur_lod];
							voxel_map[tv] = p;
						}
					}
				}
			}
		}

		// now we allocate the space and store the data in a buffer
		float *data = new float[9*element_num];
		for (map<Voxel *, std::pair<uint, uint>>::iterator it=voxel_map.begin();
				it!=voxel_map.end(); ++it){
			if(it->first->size[ctx.cur_lod]>0){
				memcpy(data+it->second.first*9, it->first->data[ctx.cur_lod], it->first->size[ctx.cur_lod]*9*sizeof(float));
			}
		}
		// organize the data for computing
		uint *offset_size = new uint[4*pair_num];
		int index = 0;
		for(candidate_entry c:candidates){
			for(candidate_info &info:c.candidates){
				for(voxel_pair &vp:info.voxel_pairs){
					offset_size[4*index] = voxel_map[vp.v1].first;
					offset_size[4*index+1] = voxel_map[vp.v1].second;
					offset_size[4*index+2] = voxel_map[vp.v2].first;
					offset_size[4*index+3] = voxel_map[vp.v2].second;
					index++;
				}
			}
		}
		assert(index==pair_num);
		ctx.packing_time += logt("organizing data with %ld elements and %ld element pairs", start, element_num, element_pair_num);


		geometry_param gp;
		gp.data = data;
		gp.pair_num = pair_num;
		gp.offset_size = offset_size;
		gp.results = ctx.results;
		gp.data_size = element_num;
		computer->get_intersect(gp);

		delete []data;
		delete []offset_size;
		voxel_map.clear();
		ctx.computation_time += logt("computation for checking intersection", start);

	}
}


//utility function to calculate the distances between voxel pairs in batch
void SpatialJoin::calculate_distance(vector<candidate_entry> &candidates, query_context &ctx){
	struct timeval start = hispeed::get_cur_time();

	const int pair_num = get_pair_num(candidates);
	ctx.results = new result_container[pair_num];
	for(int i=0;i<pair_num;i++){
		ctx.results[i].result.distance = 0;
	}

	decode_data(candidates, ctx);
	ctx.decode_time += logt("decode data", start);

	if(ctx.use_aabb){
		// build the AABB tree
		for(candidate_entry &c:candidates){
			for(candidate_info &info:c.candidates){
				if(global_ctx.etype==DT_Segment){
					info.mesh_wrapper->mesh->get_aabb_tree_segment();
				}else{
					info.mesh_wrapper->mesh->get_aabb_tree_triangle();
				}
			}
			if(global_ctx.etype==DT_Segment){
				c.mesh_wrapper->mesh->get_aabb_tree_segment();
			}else{
				c.mesh_wrapper->mesh->get_aabb_tree_triangle();
			}
		}
		ctx.packing_time += logt("building aabb tree", start);

		int index = 0;
		for(candidate_entry &c:candidates){
			for(candidate_info &info:c.candidates){
				assert(info.voxel_pairs.size()==1);
				ctx.results[index++].result.distance = c.mesh_wrapper->mesh->distance_tree(info.mesh_wrapper->mesh);
			}// end for distance_candiate list
		}// end for candidates
		// clear the trees for current LOD
		for(candidate_entry &c:candidates){
			for(candidate_info &info:c.candidates){
				info.mesh_wrapper->mesh->clear_aabb_tree();
			}
			c.mesh_wrapper->mesh->clear_aabb_tree();
		}
		ctx.computation_time += logt("computation for distance computation", start);

	}else{
		fill_voxels(candidates, ctx, global_ctx.etype);

		map<Voxel *, std::pair<uint, uint>> voxel_map;
		size_t element_num = 0;
		size_t element_pair_num = 0;

		for(candidate_entry &c:candidates){
			HiMesh_Wrapper *wrapper1 = c.mesh_wrapper;
			for(candidate_info &info:c.candidates){
				for(voxel_pair &vp:info.voxel_pairs){
					assert(vp.v1&&vp.v2);
					// update the voxel map
					for(int i=0;i<2;i++){
						Voxel *tv = i==0?vp.v1:vp.v2;
						element_pair_num += vp.v1->size[ctx.cur_lod]*vp.v2->size[ctx.cur_lod];

						if(voxel_map.find(tv)==voxel_map.end()){
							std::pair<uint, uint> p;
							p.first = element_num;
							p.second = tv->size[ctx.cur_lod];
							element_num += tv->size[ctx.cur_lod];
							voxel_map[tv] = p;
						}
					}
				}
			}
		}
		assert(element_num>0 && "there should be elements in voxels need be calculated");

		const int element_size = global_ctx.etype==DT_Segment?6:9;
		// now we allocate the space and store the data in a buffer
		float *data = new float[element_size*element_num];
		uint *offset_size = new uint[4*pair_num];

		for (map<Voxel *, std::pair<uint, uint>>::iterator it=voxel_map.begin();
				it!=voxel_map.end(); ++it){
			std::memcpy((void *)(data+it->second.first*element_size),
					(const void *)(it->first->data[ctx.cur_lod]),
					(size_t)it->first->size[ctx.cur_lod]*element_size*sizeof(float));
		}
		// organize the data for computing
		int index = 0;
		for(candidate_entry &c:candidates){
			for(candidate_info &info:c.candidates){
				for(voxel_pair &vp:info.voxel_pairs){
					assert(vp.v1!=vp.v2);
					offset_size[4*index] = voxel_map[vp.v1].first;
					offset_size[4*index+1] = voxel_map[vp.v1].second;
					offset_size[4*index+2] = voxel_map[vp.v2].first;
					offset_size[4*index+3] = voxel_map[vp.v2].second;
					index++;
				}
			}
		}
		assert(index==pair_num);
		ctx.packing_time += logt("organizing data with %ld elements and %ld element pairs", start, element_num, element_pair_num);

		geometry_param gp;
		gp.data = data;
		gp.pair_num = pair_num;
		gp.offset_size = offset_size;
		gp.results = ctx.results;
		gp.data_size = element_num;
		computer->get_distance(gp);
		delete []data;
		delete []offset_size;
		ctx.computation_time += logt("computation for distance computation", start);
	}
}

class nn_param{
public:
	pthread_mutex_t lock;
	queue<pair<Tile *, Tile *>> tile_queue;
	SpatialJoin *joiner;
	query_context ctx;
	nn_param(){
		pthread_mutex_init (&lock, NULL);
		joiner = NULL;
	}
};

void *join_unit(void *param){
	struct nn_param *nnparam = (struct nn_param *)param;
	while(!nnparam->tile_queue.empty()){
		pthread_mutex_lock(&nnparam->lock);
		if(nnparam->tile_queue.empty()){
			pthread_mutex_unlock(&nnparam->lock);
			break;
		}
		pair<Tile *, Tile *> p = nnparam->tile_queue.front();
		nnparam->tile_queue.pop();
		pthread_mutex_unlock(&nnparam->lock);
		nnparam->ctx.tile1 = p.first;
		nnparam->ctx.tile2 = p.second;

		// can only be in one of those three queries for now
		if(nnparam->ctx.query_type=="intersect"){
			nnparam->joiner->intersect(nnparam->ctx);
		}else if(nnparam->ctx.query_type=="nn"){
			nnparam->joiner->nearest_neighbor(nnparam->ctx);
		}else{
			nnparam->joiner->within(nnparam->ctx);
		}
		if(p.second!=p.first){
			delete p.second;
		}
		delete p.first;
		log("%d tile pairs left for processing",nnparam->tile_queue.size());
	}
	return NULL;
}

void SpatialJoin::join(vector<pair<Tile *, Tile *>> &tile_pairs){
	struct timeval start = hispeed::get_cur_time();
	struct nn_param param;
	for(pair<Tile *, Tile *> &p:tile_pairs){
		param.tile_queue.push(p);
	}
	param.joiner = this;
	param.ctx = global_ctx;
	pthread_t threads[global_ctx.num_thread];
	for(int i=0;i<global_ctx.num_thread;i++){
		pthread_create(&threads[i], NULL, join_unit, (void *)&param);
	}
	while(!param.tile_queue.empty()){
		usleep(10);
	}
	for(int i = 0; i < global_ctx.num_thread; i++){
		void *status;
		pthread_join(threads[i], &status);
	}
	global_ctx.report(get_time_elapsed(start));
}

}


