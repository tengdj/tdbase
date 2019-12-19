/*
 * SpatialJoin.cpp
 *
 *  Created on: Nov 11, 2019
 *      Author: teng
 */

#include <math.h>
#include <map>
#include <tuple>
#include "SpatialJoin.h"

using namespace std;

namespace hispeed{


inline bool update_voxel_pair_list(vector<voxel_pair> &voxel_pairs, range &d){
	int voxel_pair_size = voxel_pairs.size();
	for(int j=0;j<voxel_pair_size;){
		if(d>voxel_pairs[j].dist){
			return false;
		}else if(d<voxel_pairs[j].dist){
			// evict this voxel pairs
			voxel_pairs.erase(voxel_pairs.begin()+j);
			voxel_pair_size--;
		}else{
			j++;
		}
	}
	return true;
}

inline bool update_candidate_list(vector<candidate_info> &candidate_list, range &d){
	int list_size = candidate_list.size();
	for(int i=0;i<list_size;){
		if(d>candidate_list[i].distance){
			// should not keep since there is a closer one
			// in the candidate list
			return false;
		}else if(d<candidate_list[i].distance){
			// one candidate in the list cannot be the closest
			// remove it, together with the voxels, from the candidate list
			candidate_list[i].voxel_pairs.clear();
			candidate_list.erase(candidate_list.begin()+i);
			list_size--;
		}else{
			// go deep into checking the voxel pairs
			if(!update_voxel_pair_list(candidate_list[i].voxel_pairs, d)){
				return false;
			}
			// current candidate can be removed after comparing the voxel pairs
			// if all voxel pairs are farther compared to current one
			if(candidate_list[i].voxel_pairs.size()==0){
				candidate_list.erase(candidate_list.begin()+i);
				list_size--;
			}else{
				i++;
			}
		}//end if
	}//end for
	// the target one should be kept
	return true;
}

inline size_t get_pair_num(vector<candidate_entry> &candidates){
	size_t pair_num = 0;
	for(candidate_entry p:candidates){
		for(candidate_info c:p.second){
			pair_num += c.voxel_pairs.size();
		}
	}
	return pair_num;
}

bool compare_pair(pair<int, range> a1, pair<int, range> a2){
	return a1.first<a2.first;
}

void report_candidate(vector<candidate_entry> &candidates){
	for(vector<candidate_entry>::iterator it=candidates.begin();it!=candidates.end();){
		if(it->second.size()<2){
			it = candidates.erase(it);
		}else{
			it++;
		}
	}

}

vector<candidate_entry> SpatialJoin::mbb_distance(Tile *tile1, Tile *tile2){
	vector<candidate_entry> candidates;
	vector<pair<int, range>> candidate_ids;
	OctreeNode *tree = tile2->build_octree(400);
	for(int i=0;i<tile1->num_objects();i++){
		vector<candidate_info> candidate_list;
		HiMesh_Wrapper *wrapper1 = tile1->get_mesh_wrapper(i);
		tree->query_distance(&(wrapper1->box), candidate_ids);
		if(candidate_ids.empty()){
			continue;
		}
		std::sort(candidate_ids.begin(), candidate_ids.end(), compare_pair);
		int former = -1;
		for(pair<int, range> &p:candidate_ids){
			if(p.first==former){
				// duplicate
				continue;
			}
			HiMesh_Wrapper *wrapper2 = tile2->get_mesh_wrapper(p.first);
			// we firstly use the distance between the mbbs
			// of those two objects as a filter to see if this one is
			// a suitable candidate, and then we further go
			// through the voxels in two objects to shrink
			// the candidate list in a fine grained
			if(update_candidate_list(candidate_list, p.second)){
				candidate_info ci;
				for(Voxel *v1:wrapper1->voxels){
					for(Voxel *v2:wrapper2->voxels){
						range tmpd = v1->box.distance(v2->box);
						// no voxel pair in the lists is nearer
						if(update_voxel_pair_list(ci.voxel_pairs, tmpd) &&
						   update_candidate_list(candidate_list, tmpd)){
							ci.voxel_pairs.push_back(voxel_pair(v1, v2, tmpd));
							ci.distance.update(tmpd);
						}
					}
				}
				// some voxel pairs need be further evaluated
				if(ci.voxel_pairs.size()>0){
					ci.distance = p.second;
					ci.mesh_wrapper = wrapper2;
					candidate_list.push_back(ci);
				}
			}
			former = p.first;
		}
		// save the candidate list
		candidates.push_back(candidate_entry(wrapper1, candidate_list));
		candidate_ids.clear();
	}
	return candidates;
}

/*
 * the main function for getting the nearest neighbor
 *
 * */
void SpatialJoin::nearest_neighbor(Tile *tile1, Tile *tile2){
	struct timeval start = get_cur_time();

	// filtering with MBBs to get the candidate list
	vector<candidate_entry> candidates = mbb_distance(tile1, tile2);
	logt("comparing mbbs", start);
	report_candidate(candidates);

	// now we start to get the distances with progressive level of details
	int lod = base_lod;
	while(true){
		struct timeval iter_start = get_cur_time();
		const int pair_num = get_pair_num(candidates);
		if(pair_num==0){
			break;
		}
		uint *offset_size = new uint[4*pair_num];
		float *distances = new float[pair_num];
		log("%ld polyhedron has %d candidates", candidates.size(), pair_num);
		// retrieve the necessary meshes
		map<Voxel *, std::pair<uint, uint>> voxel_map;
		uint segment_num = 0;
		size_t segment_pair_num = 0;
		for(candidate_entry c:candidates){
			// the nearest neighbor is found
			if(c.second.size()<=1){
				continue;
			}
			HiMesh_Wrapper *wrapper1 = c.first;
			for(candidate_info info:c.second){
				HiMesh_Wrapper *wrapper2 = info.mesh_wrapper;
				for(voxel_pair vp:info.voxel_pairs){
					assert(vp.v1&&vp.v2);
					// not filled yet
					if(vp.v1->data.find(lod)==vp.v1->data.end()){
						// ensure the msh is extracted
						tile1->decode_to(wrapper1->id, lod);
						wrapper1->fill_voxels(DT_Segment);
					}
					if(vp.v2->data.find(lod)==vp.v2->data.end()){
						tile2->decode_to(wrapper2->id, lod);
						wrapper2->fill_voxels(DT_Segment);
					}

					// update the voxel map
					for(int i=0;i<2;i++){
						Voxel *tv = i==0?vp.v1:vp.v2;
						if(voxel_map.find(tv)==voxel_map.end()){
							std::pair<uint, uint> p;
							p.first = segment_num;
							p.second = tv->size[lod];
							segment_num += tv->size[lod];
							voxel_map[tv] = p;
						}
					}
					segment_pair_num += vp.v1->size[lod]*vp.v2->size[lod];
				}// end for voxel_pairs
			}// end for distance_candiate list
		}// end for candidates
		logt("decoded %ld voxels with %d segments %ld segment pairs for lod %d",
				start, voxel_map.size(), segment_num, segment_pair_num, lod);
		if(segment_pair_num==0){
			log("no segments is filled in this round");
			voxel_map.clear();
			continue;
		}

		// now we allocate the space and store the data in a buffer
		float *data = new float[6*segment_num];
		for (map<Voxel *, std::pair<uint, uint>>::iterator it=voxel_map.begin();
				it!=voxel_map.end(); ++it){
			memcpy(data+it->second.first*6, it->first->data[lod], it->first->size[lod]*6*sizeof(float));
		}
		// organize the data for computing
		int index = 0;
		for(candidate_entry c:candidates){
			for(candidate_info info:c.second){
				for(voxel_pair vp:info.voxel_pairs){
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
		logt("organizing data", start);
		geometry_param gp;
		gp.data = data;
		gp.pair_num = pair_num;
		gp.offset_size = offset_size;
		gp.distances = distances;
		gp.data_size = segment_num;
		computer->get_distance(gp);
		logt("get distance", start);

		// now update the distance range with the new distances
		index = 0;
		for(candidate_entry &ce:candidates){
			range min_candidate;
			min_candidate.farthest = DBL_MAX;
			for(candidate_info &ci:ce.second){
				for(voxel_pair &vp:ci.voxel_pairs){
					// update the distance
					if(vp.v1->size[lod]>0&&vp.v2->size[lod]>0){
						range dist = vp.dist;
						if(lod==100){
							// now we have a precise distance
							dist.closest = distances[index];
							dist.farthest = distances[index];
						}else{
							dist.farthest = std::min(dist.farthest, distances[index]);
						}
						vp.dist = dist;
						if(min_candidate.farthest>dist.farthest){
							min_candidate = dist;
						}
					}
					index++;
				}
			}
			update_candidate_list(ce.second, min_candidate);
		}
		report_candidate(candidates);
		logt("update candidate list", start);

		delete data;
		delete offset_size;
		delete distances;
		voxel_map.clear();
		logt("current iteration", iter_start);
		if(lod==100){
			break;
		}else{
			lod+=lod_gap;
			if(lod>100){
				lod = 100;
			}
		}
	}
}

/*
 *
 * for doing intersection
 *
 * */

inline void update_candidate_list_intersect(vector<candidate_entry> &candidates){
	for(vector<candidate_entry>::iterator it = candidates.begin();it!=candidates.end();){
		bool intersected = false;
		for(candidate_info &info:it->second){
			for(voxel_pair &vp:info.voxel_pairs){
				// if any voxel pair is ensured to be intersected
				if(vp.intersect){
					intersected = true;
					break;
				}
			}
			if(intersected){
				break;
			}
		}
		if(intersected){
			for(candidate_info &info:it->second){
				info.voxel_pairs.clear();
			}
			it->second.clear();
			it = candidates.erase(it);
		}else{
			it++;
		}
	}
}

vector<candidate_entry> SpatialJoin::mbb_intersect(Tile *tile1, Tile *tile2){
	vector<candidate_entry> candidates;
	OctreeNode *tree = tile2->build_octree(400);
	vector<int> candidate_ids;
	for(int i=0;i<tile1->num_objects();i++){
		vector<candidate_info> candidate_list;
		HiMesh_Wrapper *wrapper1 = tile1->get_mesh_wrapper(i);
		tree->query_intersect(&(wrapper1->box), candidate_ids);
		if(candidate_ids.empty()){
			continue;
		}
		std::sort(candidate_ids.begin(), candidate_ids.end());
		int former = -1;
		for(int tile2_id:candidate_ids){
			if(tile2_id==former){
				// duplicate
				continue;
			}
			HiMesh_Wrapper *wrapper2 = tile2->get_mesh_wrapper(tile2_id);
			candidate_info ci;
			ci.mesh_wrapper = wrapper2;
			for(Voxel *v1:wrapper1->voxels){
				for(Voxel *v2:wrapper2->voxels){
					if(v1->box.intersect(v2->box)){
						// a candidate not sure
						ci.voxel_pairs.push_back(voxel_pair(v1, v2));
					}
				}
			}
			// some voxel pairs need be further evaluated
			if(ci.voxel_pairs.size()>0){
				candidate_list.push_back(ci);
			}
			former = tile2_id;
		}
		candidate_ids.clear();
		// save the candidate list
		candidates.push_back(candidate_entry(wrapper1, candidate_list));
	}
	delete tree;
	candidate_ids.clear();

	return candidates;
}

/*
 * the main function for detect the intersection
 * relationship among polyhedra in the tile
 *
 * */
void SpatialJoin::intersect(Tile *tile1, Tile *tile2){
	struct timeval start = get_cur_time();

	// filtering with MBBs to get the candidate list
	vector<candidate_entry> candidates = mbb_intersect(tile1, tile2);
	logt("comparing mbbs", start);
	// evaluate the candidate list, report and remove the results confirmed
	update_candidate_list_intersect(candidates);
	logt("update candidate list", start);

	// now we start to ensure the intersection with progressive level of details
	int lod = base_lod;
	while(true){
		struct timeval iter_start = start;
		size_t pair_num = get_pair_num(candidates);
		if(pair_num==0){
			break;
		}
		log("%ld polyhedron has %ld candidates", candidates.size(), pair_num);
		// retrieve the necessary meshes
		map<Voxel *, std::pair<uint, uint>> voxel_map;
		size_t triangle_num = 0;
		for(candidate_entry c:candidates){
			HiMesh_Wrapper *wrapper1 = c.first;
			for(candidate_info info:c.second){
				HiMesh_Wrapper *wrapper2 = info.mesh_wrapper;
				for(voxel_pair vp:info.voxel_pairs){
					// not filled yet
					if(vp.v1->data.find(lod)==vp.v1->data.end()){
						tile1->decode_to(wrapper1->id, lod);
						wrapper1->fill_voxels(DT_Triangle);
					}
					if(vp.v2->data.find(lod)==vp.v2->data.end()){
						tile2->decode_to(wrapper2->id, lod);
						wrapper2->fill_voxels(DT_Triangle);
					}

					// update the voxel map
					for(int i=0;i<2;i++){
						Voxel *tv = i==0?vp.v1:vp.v2;
						if(voxel_map.find(tv)==voxel_map.end()){
							std::pair<uint, uint> p;
							p.first = triangle_num;
							p.second = tv->size[lod];
							triangle_num += tv->size[lod];
							voxel_map[tv] = p;
						}
					}
				}// end for voxel_pairs
			}// end for distance_candiate list
		}// end for candidates
		logt("decoded %ld voxels with %ld triangles for lod %d", start, voxel_map.size(), triangle_num, lod);

		// now we allocate the space and store the data in a buffer
		float *data = new float[9*triangle_num];
		for (map<Voxel *, std::pair<uint, uint>>::iterator it=voxel_map.begin();
				it!=voxel_map.end(); ++it){
			if(it->first->size[lod]>0){
				memcpy(data+it->second.first*9, it->first->data[lod], it->first->size[lod]*9*sizeof(float));
			}
		}
		// organize the data for computing
		uint *offset_size = new uint[4*pair_num];
		bool *intersect_status = new bool[pair_num];
		int index = 0;
		for(candidate_entry c:candidates){
			for(candidate_info info:c.second){
				for(voxel_pair vp:info.voxel_pairs){
					offset_size[4*index] = voxel_map[vp.v1].first;
					offset_size[4*index+1] = voxel_map[vp.v1].second;
					offset_size[4*index+2] = voxel_map[vp.v2].first;
					offset_size[4*index+3] = voxel_map[vp.v2].second;
					index++;
				}
			}
		}
		assert(index==pair_num);
		logt("organizing data", start);
		geometry_param gp;
		gp.data = data;
		gp.pair_num = pair_num;
		gp.offset_size = offset_size;
		gp.intersect = intersect_status;
		computer->get_intersect(gp);
		logt("checking intersection", start);

		// now update the intersection status and update the all candidate list
		// report results if necessary
		index = 0;
		for(int i=0;i<candidates.size();i++){
			for(int j=0;j<candidates[i].second.size();j++){
				for(int t=0;t<candidates[i].second[j].voxel_pairs.size();t++){
					// update the status
					candidates[i].second[j].voxel_pairs[t].intersect |= intersect_status[index++];
				}
			}
		}
		update_candidate_list_intersect(candidates);
		logt("update candidate list", start);

		delete data;
		delete offset_size;
		delete intersect_status;
		voxel_map.clear();

		logt("current iteration", iter_start);
		if(lod==100){
			break;
		}else{
			lod+=lod_gap;
			if(lod>100){
				lod = 100;
			}
		}
	}
}


class nn_param{
public:
	pthread_mutex_t lock;
	queue<pair<Tile *, Tile *>> tile_queue;
	SpatialJoin *joiner;
	nn_param(){
		pthread_mutex_init (&lock, NULL);
		joiner = NULL;
	}
};

void *nearest_neighbor_single(void *param){
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
		nnparam->joiner->nearest_neighbor(p.first, p.second);
	}
	return NULL;
}

void SpatialJoin::nearest_neighbor_batch(vector<pair<Tile *, Tile *>> &tile_pairs, int num_threads){
	struct nn_param param;
	for(pair<Tile *, Tile *> &p:tile_pairs){
		param.tile_queue.push(p);
	}
	param.joiner = this;
	pthread_t threads[num_threads];
	for(int i=0;i<num_threads;i++){
		pthread_create(&threads[i], NULL, nearest_neighbor_single, (void *)&param);
	}
	while(!param.tile_queue.empty()){
		usleep(10);
	}
	for(int i = 0; i < num_threads; i++){
		void *status;
		pthread_join(threads[i], &status);
	}
}

void *intersect_single(void *param){
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
		nnparam->joiner->intersect(p.first, p.second);
	}
	return NULL;
}

void SpatialJoin::intersect_batch(vector<pair<Tile *, Tile *>> &tile_pairs, int num_threads){
	struct nn_param param;
	for(pair<Tile *, Tile *> &p:tile_pairs){
		param.tile_queue.push(p);
	}
	param.joiner = this;
	pthread_t threads[num_threads];
	for(int i=0;i<num_threads;i++){
		pthread_create(&threads[i], NULL, intersect_single, (void *)&param);
	}
	for(int i = 0; i < num_threads; i++){
		void *status;
		pthread_join(threads[i], &status);
	}
}


}


