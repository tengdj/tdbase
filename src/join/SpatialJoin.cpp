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

typedef struct candidate_info_{
	HiMesh_Wrapper *mesh_wrapper;
	range distance;
	vector<std::tuple<Voxel*, Voxel *, range>> voxel_pairs;
}candidate_info;

inline bool update_voxel_pair_list(vector<std::tuple<Voxel*, Voxel *, range>> &voxel_pairs, range &d){
	int voxel_pair_size = voxel_pairs.size();
	for(int j=0;j<voxel_pair_size;){
		range cur_d = get<2>(voxel_pairs[j]);
		if(d>cur_d){
			return false;
		}else if(d<cur_d){
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

inline int get_pair_num(vector<std::pair<HiMesh_Wrapper *, vector<candidate_info>>> &candidates){
	int pair_num = 0;
	for(std::pair<HiMesh_Wrapper *, vector<candidate_info>> p:candidates){
		for(candidate_info c:p.second){
			pair_num += c.voxel_pairs.size();
		}
	}
	return pair_num;
}


void SpatialJoin::formalize_computing(){
	struct timeval start = get_cur_time();

	// filtering with MBBs to get the candidate list
	vector<std::pair<HiMesh_Wrapper *, vector<candidate_info>>> candidates;
	for(int i=0;i<tile1->num_objects();i++){
		vector<candidate_info> candidate_list;
		HiMesh_Wrapper *wrapper1 = tile1->get_mesh_wrapper(i);
		for(int j=0;j<tile2->num_objects();j++){
			// avoid self comparing
			if(tile1==tile2&&i==j){
				continue;
			}
			HiMesh_Wrapper *wrapper2 = tile2->get_mesh_wrapper(j);
			// we firstly use the distance between the mbbs
			// of those two objects as a filter to see if this one is
			// a suitable candidate, and then we further go
			// through the voxels in two objects to shrink
			// the candidate list in a fine grained
			range d = wrapper1->box.distance(wrapper2->box);
			if(update_candidate_list(candidate_list, d)){
				candidate_info ci;
				ci.distance = d;
				ci.mesh_wrapper = wrapper2;
				for(Voxel *v1:wrapper1->voxels){
					for(Voxel *v2:wrapper2->voxels){
						range tmpd = v1->box.distance(v2->box);
						// no voxel pair in the lists is nearer
						if(update_voxel_pair_list(ci.voxel_pairs, tmpd) &&
						   update_candidate_list(candidate_list, tmpd)){
							ci.voxel_pairs.push_back(std::make_tuple(v1, v2, tmpd));
						}
					}
				}
				// some voxel pairs need be further evaluated
				if(ci.voxel_pairs.size()>0){
					candidate_list.push_back(ci);
				}
			}
		}
		// save the candidate list
		candidates.push_back(std::pair<HiMesh_Wrapper *, vector<candidate_info>>(wrapper1, candidate_list));
	}
	report_time("comparing mbbs", start);

	// now we start to get the distances with progressive level of details
	for(int lod=0;lod<=100;lod+=50){
		if(true){
			int candidate_list_size = candidates.size();
			for(int i=0;i<candidate_list_size;){
				// find the nearest, report and remove
				// this entry from the candidate list
				if(candidates[i].second.size()==1){
					for(candidate_info info:candidates[i].second){
						info.voxel_pairs.clear();
					}
					candidates[i].second.clear();
					candidates.erase(candidates.begin()+i);
					candidate_list_size--;
				}else{
					i++;
				}
			}
			if(candidates.size()==0){
				return;
			}
		}

		int pair_num = get_pair_num(candidates);
		struct timeval iter_start = get_cur_time();
		printf("\n%ld polyhedron has %d candidates\n", candidates.size(), pair_num);
		// retrieve the necessary meshes
		map<Voxel *, std::pair<uint, uint>> voxel_map;
		uint segment_num = 0;
		for(std::pair<HiMesh_Wrapper *, vector<candidate_info>> c:candidates){
			HiMesh_Wrapper *wrapper1 = c.first;
			for(candidate_info info:c.second){
				HiMesh_Wrapper *wrapper2 = info.mesh_wrapper;
				for(std::tuple<Voxel *, Voxel *, range> vp:info.voxel_pairs){
					Voxel *v1 = get<0>(vp);
					Voxel *v2 = get<1>(vp);
					// not filled yet
					if(v1->data==NULL){
						tile1->get_mesh(wrapper1->id, lod);
						wrapper1->fill_voxels(lod);
					}
					if(v2->data==NULL){
						tile2->get_mesh(wrapper2->id, lod);
						wrapper2->fill_voxels(lod);
					}

					// update the voxel map
					for(int i=0;i<2;i++){
						Voxel *tv = i==0?v1:v2;
						if(voxel_map.find(tv)==voxel_map.end()){
							std::pair<uint, uint> p;
							p.first = segment_num;
							p.second = tv->size;
							segment_num += tv->size;
							voxel_map[tv] = p;
						}
					}
				}// end for voxel_pairs
			}// end for candidate_info list
		}// end for candidates
		printf("decoded %ld voxels with %d segments for lod %d\n", voxel_map.size(), segment_num, lod);
		report_time("getting data for voxels", start);

		// now we allocate the space and store the data in a buffer
		float *data = new float[6*segment_num];
		for (map<Voxel *, std::pair<uint, uint>>::iterator it=voxel_map.begin();
				it!=voxel_map.end(); ++it){
			memcpy(data+it->second.first*6, it->first->data, it->first->size*6*sizeof(float));
		}
		// organize the data for computing
		uint *offset_size = new uint[4*pair_num];
		float *distances = new float[pair_num];
		int index = 0;
		for(std::pair<HiMesh_Wrapper *, vector<candidate_info>> c:candidates){
			for(candidate_info info:c.second){
				for(std::tuple<Voxel *, Voxel*, range> vp:info.voxel_pairs){
					Voxel *v1 = get<0>(vp);
					Voxel *v2 = get<1>(vp);
					assert(v1!=v2);
					offset_size[4*index] = voxel_map[v1].first;
					offset_size[4*index+1] = voxel_map[v1].second;
					offset_size[4*index+2] = voxel_map[v2].first;
					offset_size[4*index+3] = voxel_map[v2].second;
					index++;
				}
			}
		}
		assert(index==pair_num);
		report_time("organizing data", start);
		if(true){
			hispeed::SegDist_batch_gpu(data, offset_size, distances, pair_num, segment_num);
			report_time("get distance with GPU", start);
		}else{
			hispeed::SegDist_batch(data, offset_size, distances, pair_num, hispeed::get_num_threads());
			report_time("get distance with CPU", start);
		}

		// now update the distance range with the new distances
		index = 0;
		std::vector<range> tmp_distances;
		for(int i=0;i<candidates.size();i++){
			for(int j=0;j<candidates[i].second.size();j++){
				for(int t=0;t<candidates[i].second[j].voxel_pairs.size();t++){
					// update the distance
					if(get<0>(candidates[i].second[j].voxel_pairs[t])->size>0&&
					   get<1>(candidates[i].second[j].voxel_pairs[t])->size>0){
						range dist = get<2>(candidates[i].second[j].voxel_pairs[t]);
						if(lod==100){
							// now we have a precise distance
							dist.closest = distances[index];
							dist.farthest = distances[index];
						}else{
							dist.farthest = std::min(dist.farthest, distances[index]);
						}
						get<2>(candidates[i].second[j].voxel_pairs[t]) = dist;
						tmp_distances.push_back(dist);
					}
					index++;
				}
			}
			for(range r:tmp_distances){
				update_candidate_list(candidates[i].second, r);
			}
			tmp_distances.clear();
		}
		report_time("update distance range", start);

		// reset the voxels
		for(int i=0;i<candidates.size();i++){
			candidates[i].first->reset();
			for(int j=0;j<candidates[i].second.size();j++){
				candidates[i].second[j].mesh_wrapper->reset();
			}
		}

		delete data;
		delete offset_size;
		delete distances;
		voxel_map.clear();
		report_time("current iteration", iter_start, false);
		pair_num = get_pair_num(candidates);
		if(pair_num==0){
			break;
		}
	}
}


}


