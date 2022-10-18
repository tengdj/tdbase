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

//utility function to calculate the distances between voxel pairs in batch
float *SpatialJoin::calculate_distance(vector<candidate_entry> &candidates, query_context &ctx, const int lod){
	const int pair_num = get_pair_num(candidates);
	const int candidate_num = get_candidate_num(candidates);
	float *distances = new float[pair_num];
	for(int i=0;i<pair_num;i++){
		distances[i] = 0;
	}

	struct timeval start = hispeed::get_cur_time();
	if(ctx.use_aabb){
		assert(pair_num==candidate_num && "no shape-aware indexing should be applied in aabb");
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
		ctx.packing_time += hispeed::get_time_elapsed(start, false);
		logt("building aabb tree", start);

		int index = 0;
		for(candidate_entry &c:candidates){
			for(candidate_info &info:c.candidates){
				assert(info.voxel_pairs.size()==1);
				distances[index++] = c.mesh_wrapper->mesh->distance_tree(info.mesh_wrapper->mesh);
			}// end for distance_candiate list
		}// end for candidates
		// clear the trees for current LOD
		for(candidate_entry &c:candidates){
			for(candidate_info &info:c.candidates){
				info.mesh_wrapper->mesh->clear_aabb_tree();
			}
			c.mesh_wrapper->mesh->clear_aabb_tree();
		}
	}else{

		map<Voxel *, std::pair<uint, uint>> voxel_map;
		uint element_num = 0;

		for(candidate_entry &c:candidates){
			HiMesh_Wrapper *wrapper1 = c.mesh_wrapper;
			for(candidate_info &info:c.candidates){
				for(voxel_pair &vp:info.voxel_pairs){
					assert(vp.v1&&vp.v2);
					// update the voxel map
					for(int i=0;i<2;i++){
						Voxel *tv = i==0?vp.v1:vp.v2;
						if(voxel_map.find(tv)==voxel_map.end()){
							std::pair<uint, uint> p;
							p.first = element_num;
							p.second = tv->size[lod];
							element_num += tv->size[lod];
							voxel_map[tv] = p;
						}
					}
				}
			}
		}

		assert(element_num>0 && "there should be elements in voxels need be calculated");

		int element_size = global_ctx.etype==DT_Segment?6:9;
		// now we allocate the space and store the data in a buffer
		float *data = new float[element_size*element_num];
		uint *offset_size = new uint[4*pair_num];

		for (map<Voxel *, std::pair<uint, uint>>::iterator it=voxel_map.begin();
				it!=voxel_map.end(); ++it){
			std::memcpy((void *)(data+it->second.first*element_size),
					(const void *)(it->first->data[lod]),
					(size_t)it->first->size[lod]*element_size*sizeof(float));
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
		ctx.packing_time += hispeed::get_time_elapsed(start, false);
		logt("organizing data", start);
		geometry_param gp;
		gp.data = data;
		gp.pair_num = pair_num;
		gp.offset_size = offset_size;
		gp.distances = distances;
		gp.data_size = element_num;
		computer->get_distance(gp);
		delete []data;
		delete []offset_size;
	}
	return distances;
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
		// can only be in one of those three queries for now
		if(nnparam->ctx.query_type=="intersect"){
			nnparam->joiner->intersect(p.first, p.second,nnparam->ctx);
		}else if(nnparam->ctx.query_type=="nn"){
			nnparam->joiner->nearest_neighbor(p.first, p.second,nnparam->ctx);
		}else{
			nnparam->joiner->within(p.first, p.second, nnparam->ctx);
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


