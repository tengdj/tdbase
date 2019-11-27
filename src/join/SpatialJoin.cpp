/*
 * SpatialJoin.cpp
 *
 *  Created on: Nov 11, 2019
 *      Author: teng
 */


#include "SpatialJoin.h"
#include <map>
using namespace std;

namespace hispeed{

typedef struct candidate_info_{
	HiMesh_Wrapper *mesh_wrapper;
	range distance;
}candidate_info;

void SpatialJoin::formalize_computing(){
	struct timeval start = get_cur_time();

	// filtering with MBBs to get the candidate list
	vector<std::pair<int, vector<candidate_info>>> candidates;
	int pair_num = 0;
	for(int i=0;i<tile1->num_objects();i++){
		vector<candidate_info> candidate_list;
		aab b1 = tile1->get_mbb(i);
		for(int j=0;j<tile2->num_objects();j++){
			// avoid self comparing
			if(tile1==tile2&&i==j){
				continue;
			}
			HiMesh_Wrapper *wrapper = tile2->get_mesh_wrapper(j);
			range d = b1.distance(wrapper->box);
			// now update the list
			bool keep = true;
			int list_size = candidate_list.size();
			for(int i=0;i<list_size;){
				// should not keep since there is a closer one
				if(d>candidate_list[i].distance){
					keep = false;
					break;
				// one in the list cannot be the closest because of this one
				}else if(d<candidate_list[i].distance){
					candidate_list.erase(candidate_list.begin()+i);
					list_size--;
				}else{
					i++;
				}
			}
			if(keep){
				candidate_info ci;
				ci.distance = d;
				ci.mesh_wrapper = wrapper;
				candidate_list.push_back(ci);
			}
		}
		pair_num += candidate_list.size();
		// save the candidate list
		candidates.push_back(std::pair<int, vector<candidate_info>>(i, candidate_list));
	}
	report_time("comparing mbbs", start);

	// now we start to get the distances with progressive level of details
	for(int lod=0;lod<=100;lod+=50){
		printf("\n%ld polyhedron has %d candidates\n", candidates.size(), pair_num);
		// retrieve the necessary meshes
		map<HiMesh *, std::pair<uint, uint>> mesh_map;
		uint segment_num = 0;
		for(std::pair<int, vector<candidate_info>> c:candidates){
			HiMesh *mesh1 = tile1->get_mesh(c.first, lod);
			if(mesh_map.find(mesh1)==mesh_map.end()){
				std::pair<uint, uint> p;
				p.first = segment_num;
				p.second = mesh1->get_segment_num();
				segment_num += p.second;
				mesh_map[mesh1] = p;
			}
			for(candidate_info info:c.second){
				HiMesh *mesh2 = tile2->get_mesh(info.mesh_wrapper->id, lod);
				if(mesh_map.find(mesh2)==mesh_map.end()){
					std::pair<uint, uint> p;
					p.first = segment_num;
					p.second = mesh2->get_segment_num();
					segment_num += p.second;
					mesh_map[mesh2] = p;
				}
			}
		}
		printf("decoded %ld meshes with %d segments for lod %d\n", mesh_map.size(), segment_num, lod);
		report_time("getting meshes", start);

		// now we allocate the space and store the data in
		// a buffer
		float *data = new float[6*segment_num];
		for (map<HiMesh *, std::pair<uint, uint>>::iterator it=mesh_map.begin();
				it!=mesh_map.end(); ++it){
			assert(it->second.second==it->first->get_segments(data+it->second.first*6));
		}
		// organize the data for computing
		uint *offset_size = new uint[4*pair_num];
		float *distances = new float[pair_num];
		int index = 0;
		for(std::pair<int, vector<candidate_info>> c:candidates){
			HiMesh *mesh1 = tile1->get_mesh(c.first, lod);
			for(candidate_info info:c.second){
				HiMesh *mesh2 = tile2->get_mesh(info.mesh_wrapper->id, lod);
				assert(mesh1!=mesh2);
				offset_size[4*index] = mesh_map[mesh1].first;
				offset_size[4*index+1] = mesh_map[mesh1].second;
				offset_size[4*index+2] = mesh_map[mesh2].first;
				offset_size[4*index+3] = mesh_map[mesh2].second;
				index++;
			}
		}
		assert(index==pair_num);
		report_time("organizing data", start);

		hispeed::SegDist_batch_gpu(data, offset_size, distances, pair_num, segment_num);
		report_time("get distance with GPU", start);

		// now update the distance range with the new distances
		index = 0;
		pair_num = 0;
		int candidate_list_size = candidates.size();
		for(int i=0;i<candidate_list_size;){
			for(int j=0;j<candidates[i].second.size();j++){
//				assert(candidates[i].second[j].distance.farthest>=distances[index]&&
//					   candidates[i].second[j].distance.closest<=distances[index]);
				// now we have a precise distance
				if(lod==100){
					candidates[i].second[j].distance.closest = distances[index];
					candidates[i].second[j].distance.farthest = distances[index];
				}else{
					candidates[i].second[j].distance.farthest =
								std::min(candidates[i].second[j].distance.farthest,distances[index]);
				}
				index++;
			}

			int candidate_num = candidates[i].second.size();
			for(int j = 0;j<candidate_num;){
				bool keep = true;
				for(int jj=j+1;jj<candidate_num;){
					// should not keep since there is a closer one
					if(candidates[i].second[j].distance>candidates[i].second[jj].distance){
						keep = false;
						break;
					// one in the list cannot be the closest because of this one
					}else if(candidates[i].second[j].distance<candidates[i].second[jj].distance){
						candidate_num--;
						candidates[i].second.erase(candidates[i].second.begin()+jj);
					}else{
						jj++;
					}
				}
				if(!keep){
					candidate_num--;
					candidates[i].second.erase(candidates[i].second.begin()+j);
				}else{
					j++;
				}
			}
			// find the nearest
			if(candidates[i].second.size()==1){
				candidates.erase(candidates.begin()+i);
				candidate_list_size--;
			}else{
				pair_num += candidates[i].second.size();
				i++;
			}
		}
		report_time("update distance range", start);
		delete data;
		delete offset_size;
		delete distances;
		mesh_map.clear();
		if(pair_num==0){
			break;
		}
	}


}


}


