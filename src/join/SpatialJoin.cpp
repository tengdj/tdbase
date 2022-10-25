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

void SpatialJoin::clear_data(vector<candidate_entry> &candidates, query_context &ctx){
	for(candidate_entry &c:candidates){
		HiMesh_Wrapper *wrapper1 = c.mesh_wrapper;
		for(candidate_info &info:c.candidates){
			for(voxel_pair &vp:info.voxel_pairs){
				assert(vp.v1&&vp.v2);
				vp.v1->reset();
				vp.v2->reset();
			}// end for voxel_pairs
		}// end for distance_candiate list
	}// end for candidates
}

void SpatialJoin::decode_data(vector<candidate_entry> &candidates, query_context &ctx){
	for(candidate_entry &c:candidates){
		HiMesh_Wrapper *wrapper1 = c.mesh_wrapper;
		for(candidate_info &info:c.candidates){
			for(voxel_pair &vp:info.voxel_pairs){
				assert(vp.v1&&vp.v2);
				// ensure the mesh is extracted and decoded
				ctx.tile1->decode_to(wrapper1->id, ctx.cur_lod);
				ctx.tile2->decode_to(info.mesh_wrapper->id, ctx.cur_lod);
			}// end for voxel_pairs
		}// end for distance_candiate list
	}// end for candidates
}

void SpatialJoin::fill_voxels(vector<candidate_entry> &candidates, query_context &ctx){
	for(candidate_entry &c:candidates){
		HiMesh_Wrapper *wrapper1 = c.mesh_wrapper;
		for(candidate_info &info:c.candidates){
			for(voxel_pair &vp:info.voxel_pairs){
				assert(vp.v1&&vp.v2);
				// not filled yet
				// ensure the mesh is extracted
				wrapper1->fill_voxels();
				info.mesh_wrapper->fill_voxels();
			}// end for voxel_pairs
		}// end for distance_candiate list
	}// end for candidates
}

geometry_param SpatialJoin::packing_data(vector<candidate_entry> &candidates, query_context &ctx){
	geometry_param gp;
	gp.pair_num = get_pair_num(candidates);
	gp.element_num = 0;
	gp.element_pair_num = 0;
	gp.results = ctx.results;
	map<Voxel *, std::pair<uint, uint>> voxel_map;

	fill_voxels(candidates, ctx);

	for(candidate_entry &c:candidates){
		HiMesh_Wrapper *wrapper1 = c.mesh_wrapper;
		for(candidate_info &info:c.candidates){
			for(voxel_pair &vp:info.voxel_pairs){
				assert(vp.v1->data&&vp.v2->data);
				gp.element_pair_num += vp.v1->data->size*vp.v2->data->size;
				// update the voxel map
				for(int i=0;i<2;i++){
					Voxel *tv = i==0?vp.v1:vp.v2;
					if(voxel_map.find(tv)==voxel_map.end()){
						std::pair<uint, uint> p;
						p.first = gp.element_num;
						p.second = tv->data->size;
						gp.element_num += tv->data->size;
						voxel_map[tv] = p;
					}
				}
			}
		}
	}

	// now we allocate the space and store the data in the buffer
	gp.data = new float[9*gp.element_num];
	for (map<Voxel *, std::pair<uint, uint>>::iterator it=voxel_map.begin();
			it!=voxel_map.end(); ++it){
		assert(it->first->data);
		if(it->first->data->size>0){
			memcpy(gp.data+it->second.first*9, it->first->data->data, it->first->data->size*9*sizeof(float));
		}
	}
	// organize the data for computing
	gp.offset_size = new uint[4*gp.pair_num];
	int index = 0;
	for(candidate_entry c:candidates){
		for(candidate_info &info:c.candidates){
			for(voxel_pair &vp:info.voxel_pairs){
				gp.offset_size[4*index] = voxel_map[vp.v1].first;
				gp.offset_size[4*index+1] = voxel_map[vp.v1].second;
				gp.offset_size[4*index+2] = voxel_map[vp.v2].first;
				gp.offset_size[4*index+3] = voxel_map[vp.v2].second;
				index++;
			}
		}
	}
	assert(index==gp.pair_num);
	voxel_map.clear();
	return gp;
}

void SpatialJoin::check_intersection(vector<candidate_entry> &candidates, query_context &ctx){
	struct timeval start = hispeed::get_cur_time();

	const int pair_num = get_pair_num(candidates);

	ctx.results = new result_container[pair_num];
	for(int i=0;i<pair_num;i++){
		ctx.results[i].result.intersected = false;
	}
	clear_data(candidates, ctx);
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
		geometry_param gp = packing_data(candidates, ctx);
		ctx.packing_time += logt("organizing data with %ld elements and %ld element pairs", start, gp.element_num, gp.element_pair_num);

		computer->get_intersect(gp);
		delete []gp.data;
		delete []gp.offset_size;
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

	clear_data(candidates, ctx);
	decode_data(candidates, ctx);
	ctx.decode_time += logt("decode data", start);

	if(ctx.use_aabb){
		// build the AABB tree
		for(candidate_entry &c:candidates){
			for(candidate_info &info:c.candidates){
				info.mesh_wrapper->mesh->get_aabb_tree_triangle();
			}
			c.mesh_wrapper->mesh->get_aabb_tree_triangle();
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
		geometry_param gp = packing_data(candidates, ctx);
		ctx.packing_time += logt("organizing data with %ld elements and %ld element pairs", start, gp.element_num, gp.element_pair_num);
		computer->get_distance(gp);
		delete []gp.data;
		delete []gp.offset_size;
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
		nnparam->ctx.tile1->init();
		if(nnparam->ctx.tile1!=nnparam->ctx.tile2){
			nnparam->ctx.tile2->init();
		}
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


