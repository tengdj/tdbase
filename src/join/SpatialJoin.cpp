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

inline size_t get_candidate_num(vector<candidate_entry> &candidates){
	size_t candidate_num = 0;
	for(candidate_entry p:candidates){
		candidate_num += p.second.size();
	}
	return candidate_num;
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

void SpatialJoin::report_time(double t){
	cout<<"total, index, decode, packing, computation, updatelist, other"<<endl;
	cout<<t<<","
		<<t*global_index_time/global_total_time<<","
		<<t*global_decode_time/global_total_time<<","
		<<t*global_packing_time/global_total_time<<","
		<<t*global_computation_time/global_total_time<<","
		<<t*global_updatelist_time/global_total_time<<","
		<<t*(global_total_time-global_decode_time-global_computation_time-global_index_time)/global_total_time<<endl;
//	cout<<"total, decode, computation, other"<<endl;
//	t /= 1000;
//	cout<<t<<","
//		<<t*global_decode_time/global_total_time<<","
//		<<t*global_computation_time/global_total_time<<","
//		<<t*(global_total_time-global_decode_time-global_computation_time)/global_total_time<<endl;
}

vector<candidate_entry> SpatialJoin::mbb_distance(Tile *tile1, Tile *tile2, double max_dist){
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
						if(tmpd.closest>max_dist){
							continue;
						}
						// no voxel pair in the list is nearer
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
void SpatialJoin::within(Tile *tile1, Tile *tile2, double max_dist){
	struct timeval start = get_cur_time();
	struct timeval very_start = get_cur_time();

	double index_time = 0;
	double decode_time = 0;
	double packing_time = 0;
	double computation_time = 0;
	double updatelist_time = 0;
	// filtering with MBBs to get the candidate list
	vector<candidate_entry> candidates = mbb_distance(tile1, tile2, max_dist);
	index_time += get_time_elapsed(start, false);
	logt("comparing mbbs", start);

	// now we start to get the distances with progressive level of details
	if(lods.size()==0){
		for(int lod = base_lod;lod<=top_lod;lod+=lod_gap){
			lods.push_back(lod);
		}
		if(lods[lods.size()-1]<top_lod){
			lods.push_back(top_lod);
		}
	}

	for(int lod:lods){
		struct timeval iter_start = get_cur_time();
		const int pair_num = get_pair_num(candidates);
		if(pair_num==0){
			break;
		}
		uint *offset_size = new uint[4*pair_num];
		float *distances = new float[pair_num];
		size_t candidate_num = get_candidate_num(candidates);
		log("%ld polyhedron has %d candidates %f voxel pairs per candidate", candidates.size(), candidate_num, (1.0*pair_num)/candidates.size());
		// retrieve the necessary meshes
		map<Voxel *, std::pair<uint, uint>> voxel_map;
		uint segment_num = 0;
		size_t segment_pair_num = 0;

		int candidate_length = 0;
		for(candidate_entry &e:candidates){
			for(candidate_info &i:e.second){
				candidate_length += i.voxel_pairs.size();
			}
		}

		for(candidate_entry &c:candidates){
			HiMesh_Wrapper *wrapper1 = c.first;
			for(candidate_info &info:c.second){
				HiMesh_Wrapper *wrapper2 = info.mesh_wrapper;
				for(voxel_pair &vp:info.voxel_pairs){
					assert(vp.v1&&vp.v2);
					// not filled yet
					if(vp.v1->data.find(lod)==vp.v1->data.end()){
						// ensure the mesh is extracted
						tile1->decode_to(wrapper1->id, lod);
						wrapper1->fill_voxels(DT_Segment, lod==lods[lods.size()-1]);
					}
					if(vp.v2->data.find(lod)==vp.v2->data.end()){
						tile2->decode_to(wrapper2->id, lod);
						wrapper2->fill_voxels(DT_Segment, lod==lods[lods.size()-1]);
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
		decode_time += hispeed::get_time_elapsed(start, false);
		logt("decoded %ld voxels with %d segments %ld segment pairs for lod %d",
				start, voxel_map.size(), segment_num, segment_pair_num, lod);
		if(segment_pair_num==0){
			log("no segments is filled in this round");
			voxel_map.clear();
			continue;
		}
		tile1->reset_time();

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
		packing_time += hispeed::get_time_elapsed(start, false);
		logt("organizing data", start);
		geometry_param gp;
		gp.data = data;
		gp.pair_num = pair_num;
		gp.offset_size = offset_size;
		gp.distances = distances;
		gp.data_size = segment_num;
		computer->get_distance(gp);
		computation_time += hispeed::get_time_elapsed(start, false);
		logt("get distance", start);

		// now update the candidate list with the new meshes
		index = 0;
		for(candidate_entry &ce:candidates){
			bool determined = false;
			for(vector<candidate_info>::iterator it=ce.second.begin();it!=ce.second.end();){
				for(voxel_pair &vp:it->voxel_pairs){
					// update the distance
					if(!determined&&vp.v1->size[lod]>0&&vp.v2->size[lod]>0){
						range dist = vp.dist;
						if(lod==lods[lods.size()-1]){
							// now we have a precise distance
							dist.closest = distances[index];
							dist.farthest = distances[index];
						}else{
							dist.farthest = std::min(dist.farthest, distances[index]);
						}
						vp.dist = dist;
						if(dist.farthest<max_dist){
							determined = true;
						}
					}
					index++;
				}
				if(determined){
					it = ce.second.erase(it);
				}else{
					it++;
				}
			}
		}
		updatelist_time += hispeed::get_time_elapsed(start, false);
		logt("update candidate list", start);

		delete data;
		delete offset_size;
		delete distances;
		voxel_map.clear();
		logt("current iteration", iter_start);

		if(lod==lods[lods.size()-1]){
			break;
		}
	}
	pthread_mutex_lock(&g_lock);
	global_index_time += index_time;
	global_decode_time += decode_time;
	global_packing_time += packing_time;
	global_computation_time += computation_time;
	global_updatelist_time += updatelist_time;
	global_total_time += hispeed::get_time_elapsed(very_start, false);
	pthread_mutex_unlock(&g_lock);

}


/*
 * the main function for getting the nearest neighbor
 *
 * */
void SpatialJoin::nearest_neighbor(Tile *tile1, Tile *tile2){
	struct timeval start = get_cur_time();
	struct timeval very_start = get_cur_time();

	double index_time = 0;
	double decode_time = 0;
	double packing_time = 0;
	double computation_time = 0;
	double updatelist_time = 0;
	// filtering with MBBs to get the candidate list
	vector<candidate_entry> candidates = mbb_distance(tile1, tile2);
	index_time += get_time_elapsed(start, false);
	logt("comparing mbbs", start);
	report_candidate(candidates);


	// now we start to get the distances with progressive level of details
	if(lods.size()==0){
		for(int lod = base_lod;lod<=top_lod;lod+=lod_gap){
			lods.push_back(lod);
		}
		if(lods[lods.size()-1]<top_lod){
			lods.push_back(top_lod);
		}
	}


	for(int lod:lods){
		struct timeval iter_start = get_cur_time();
		const int pair_num = get_pair_num(candidates);
		if(pair_num==0){
			break;
		}
		uint *offset_size = new uint[4*pair_num];
		float *distances = new float[pair_num];
		size_t candidate_num = get_candidate_num(candidates);
		log("%ld polyhedron has %d candidates %f voxel pairs per candidate", candidates.size(), candidate_num, (1.0*pair_num)/candidates.size());
		// retrieve the necessary meshes
		map<Voxel *, std::pair<uint, uint>> voxel_map;
		uint segment_num = 0;
		size_t segment_pair_num = 0;

		int candidate_length = 0;
		for(candidate_entry &e:candidates){
			for(candidate_info &i:e.second){
				candidate_length += i.voxel_pairs.size();
			}
		}

		for(candidate_entry &c:candidates){
			// the nearest neighbor is found
			if(c.second.size()<=1){
				continue;
			}
			HiMesh_Wrapper *wrapper1 = c.first;
			for(candidate_info &info:c.second){
				HiMesh_Wrapper *wrapper2 = info.mesh_wrapper;
				for(voxel_pair &vp:info.voxel_pairs){
					assert(vp.v1&&vp.v2);
					// not filled yet
					if(vp.v1->data.find(lod)==vp.v1->data.end()){
						// ensure the mesh is extracted
						tile1->decode_to(wrapper1->id, lod);
						wrapper1->fill_voxels(DT_Segment, lod==lods[lods.size()-1]);
					}
					if(vp.v2->data.find(lod)==vp.v2->data.end()){
						tile2->decode_to(wrapper2->id, lod);
						wrapper2->fill_voxels(DT_Segment, lod==lods[lods.size()-1]);
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
		decode_time += hispeed::get_time_elapsed(start, false);
		logt("decoded %ld voxels with %d segments %ld segment pairs for lod %d",
				start, voxel_map.size(), segment_num, segment_pair_num, lod);
		if(segment_pair_num==0){
			log("no segments is filled in this round");
			voxel_map.clear();
			continue;
		}
//		cerr<<"\ndecoding time\t"<<tile1->decode_time
//			<<"\n\tretrieve time\t"<< tile1->retrieve_time
//			<<"\n\t\tdisk time\t" << tile1->disk_time
//			<<"\n\t\tmalloc time\t"<<tile1->malloc_time
//			<<"\n\t\tnewmesh time\t"<<tile1->newmesh_time
//			<<"\n\tadvance time\t"<< tile1->advance_time
//			<<"\nfilling time\t"<<fill_time
//			<<endl<<endl;
		tile1->reset_time();

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
		packing_time += hispeed::get_time_elapsed(start, false);
		logt("organizing data", start);
		geometry_param gp;
		gp.data = data;
		gp.pair_num = pair_num;
		gp.offset_size = offset_size;
		gp.distances = distances;
		gp.data_size = segment_num;
		computer->get_distance(gp);
		computation_time += hispeed::get_time_elapsed(start, false);
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
						if(lod==lods[lods.size()-1]){
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
		updatelist_time += hispeed::get_time_elapsed(start, false);
		logt("update candidate list", start);

		delete data;
		delete offset_size;
		delete distances;
		voxel_map.clear();
		logt("current iteration", iter_start);

		if(lod==lods[lods.size()-1]){
			break;
		}
	}
	pthread_mutex_lock(&g_lock);
	global_index_time += index_time;
	global_decode_time += decode_time;
	global_packing_time += packing_time;
	global_computation_time += computation_time;
	global_updatelist_time += updatelist_time;
	global_total_time += hispeed::get_time_elapsed(very_start, false);
	pthread_mutex_unlock(&g_lock);

}




void SpatialJoin::nearest_neighbor_aabb(Tile *tile1, Tile *tile2){
	struct timeval start = get_cur_time();
	struct timeval very_start = get_cur_time();
	double index_time = 0;
	double decode_time = 0;
	double packing_time = 0;
	double computation_time = 0;
	// filtering with MBBs to get the candidate list
	vector<candidate_entry> candidates = mbb_distance(tile1, tile2);
	index_time += hispeed::get_time_elapsed(start, false);
	logt("comparing mbbs", start);


	// generate the aabb tree for all the referred polyhedrons
	int candidate_size = 0;
	for(candidate_entry c:candidates){
		if(c.second.size()<=1){
			continue;
		}
		for(candidate_info info:c.second){
			HiMesh_Wrapper *wrapper2 = info.mesh_wrapper;
			tile2->decode_to(wrapper2->id,100);
		}
		tile1->decode_to(c.first->id, 100);
		candidate_size += c.second.size();
	}
	decode_time += hispeed::get_time_elapsed(start, false);
	logt("decode data",start);
	log("%ld polyhedron has %d candidates", candidates.size(), candidate_size);

	for(candidate_entry c:candidates){
		if(c.second.size()<=1){
			continue;
		}
		for(candidate_info info:c.second){
			HiMesh_Wrapper *wrapper2 = info.mesh_wrapper;
			wrapper2->mesh->get_aabb_tree();
		}
	}
	packing_time += hispeed::get_time_elapsed(start, false);
	logt("build aabb tree",start);

	for(candidate_entry c:candidates){
		// the nearest neighbor is found
		if(c.second.size()<=1){
			continue;
		}
		HiMesh_Wrapper *wrapper1 = c.first;
		double min_dist = DBL_MAX;
		vector<Point> vertices;
		wrapper1->mesh->get_vertices(vertices);
		for(candidate_info info:c.second){
			HiMesh_Wrapper *wrapper2 = info.mesh_wrapper;
			int index = 0;
			for(Point &p:vertices){
				FT sqd = wrapper2->mesh->get_aabb_tree()->squared_distance(p);
				double distance = (double)CGAL::to_double(sqd);
				if(min_dist>distance){
					min_dist = distance;
				}
			}
		}// end for distance_candiate list
		vertices.clear();
	}// end for candidates
	computation_time += hispeed::get_time_elapsed(start, false);
	logt("getting distance", start);

	pthread_mutex_lock(&g_lock);
	global_index_time += index_time;
	global_decode_time += decode_time;
	global_packing_time += packing_time;
	global_computation_time += computation_time;
	global_total_time += hispeed::get_time_elapsed(very_start, false);
	pthread_mutex_unlock(&g_lock);

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
 * the main function for detecting the intersection
 * relationship among polyhedra in the tile
 *
 * */
void SpatialJoin::intersect(Tile *tile1, Tile *tile2){
	struct timeval start = get_cur_time();
	struct timeval very_start = get_cur_time();
	double index_time = 0;
	double decode_time = 0;
	double packing_time = 0;
	double updatelist_time = 0;
	double computation_time = 0;

	// filtering with MBBs to get the candidate list
	vector<candidate_entry> candidates = mbb_intersect(tile1, tile2);
	index_time += hispeed::get_time_elapsed(start,false);
	logt("comparing mbbs", start);
	// evaluate the candidate list, report and remove the results confirmed
	update_candidate_list_intersect(candidates);
	updatelist_time += hispeed::get_time_elapsed(start, false);
	logt("update candidate list", start);

	// now we start to ensure the intersection with progressive level of details
	size_t triangle_pair_num = 0;
	double fill_time = 0;
	if(lods.size()==0){
		for(int lod = base_lod;lod<=top_lod;lod+=lod_gap){
			lods.push_back(lod);
		}
		if(lods[lods.size()-1]<top_lod){
			lods.push_back(top_lod);
		}
	}
	for(int lod:lods){
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
						timeval cur = hispeed::get_cur_time();
						wrapper1->fill_voxels(DT_Triangle, lod==lods[lods.size()-1]);
						fill_time += hispeed::get_time_elapsed(cur, true);
					}
					if(vp.v2->data.find(lod)==vp.v2->data.end()){
						tile2->decode_to(wrapper2->id, lod);
						timeval cur = hispeed::get_cur_time();
						wrapper2->fill_voxels(DT_Triangle, lod==lods[lods.size()-1]);
						fill_time += hispeed::get_time_elapsed(cur, true);
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
					triangle_pair_num += vp.v1->size[lod]*vp.v2->size[lod];
				}// end for voxel_pairs
			}// end for distance_candiate list
		}// end for candidates
		decode_time += hispeed::get_time_elapsed(start, false);
		logt("decoded %ld voxels with %ld triangles %ld pairs for lod %d",
				start, voxel_map.size(), triangle_num, triangle_pair_num, lod);

//		cerr<<"\ndecoding time\t"<<tile1->decode_time+tile2->decode_time
//			<<"\n\tretrieve time\t"<< tile1->retrieve_time+tile2->retrieve_time
//			<<"\n\t\tdisk time\t" << tile1->disk_time+tile2->disk_time
//			<<"\n\t\tmalloc time\t"<<tile1->malloc_time+tile2->malloc_time
//			<<"\n\t\tnewmesh time\t"<<tile1->newmesh_time+tile2->newmesh_time
//			<<"\n\tadvance time\t"<< tile1->advance_time+tile2->advance_time
//			<<"\nfilling time\t"<<fill_time
//			<<endl<<endl;
		tile1->reset_time();
		tile2->reset_time();
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
		for(int i=0;i<pair_num;i++){
			intersect_status[i] = false;
		}
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
		packing_time += hispeed::get_time_elapsed(start, false);
		logt("organizing data", start);
		geometry_param gp;
		gp.data = data;
		gp.pair_num = pair_num;
		gp.offset_size = offset_size;
		gp.intersect = intersect_status;
		computer->get_intersect(gp);
		computation_time += hispeed::get_time_elapsed(start, false);
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
		updatelist_time += hispeed::get_time_elapsed(start, false);
		logt("update candidate list", start);

		delete []data;
		delete []offset_size;
		delete []intersect_status;
		voxel_map.clear();

		logt("current iteration", iter_start);
	}

	pthread_mutex_lock(&g_lock);
	global_index_time += index_time;
	global_decode_time += decode_time;
	global_packing_time += packing_time;
	global_computation_time += computation_time;
	global_updatelist_time += updatelist_time;
	global_total_time += hispeed::get_time_elapsed(very_start, false);
	pthread_mutex_unlock(&g_lock);

}


class nn_param{
public:
	double max_dist = DBL_MAX;
	pthread_mutex_t lock;
	queue<pair<Tile *, Tile *>> tile_queue;
	SpatialJoin *joiner;
	nn_param(){
		pthread_mutex_init (&lock, NULL);
		joiner = NULL;
	}
	bool ispeed = false;
};

void *within_single(void *param){
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
		nnparam->joiner->within(p.first, p.second, nnparam->max_dist);
		if(p.second!=p.first){
			delete p.second;
		}
		delete p.first;
	}
	return NULL;
}

void SpatialJoin::within_batch(vector<pair<Tile *, Tile *>> &tile_pairs, int num_threads, double max_dist){
	struct nn_param param;
	for(pair<Tile *, Tile *> &p:tile_pairs){
		param.tile_queue.push(p);
	}
	param.joiner = this;
	param.max_dist = max_dist;
	pthread_t threads[num_threads];
	for(int i=0;i<num_threads;i++){
		pthread_create(&threads[i], NULL, within_single, (void *)&param);
	}
	while(!param.tile_queue.empty()){
		usleep(10);
	}
	for(int i = 0; i < num_threads; i++){
		void *status;
		pthread_join(threads[i], &status);
	}
}


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
		if(nnparam->ispeed){
			p.first->disable_innerpart();
			p.second->disable_innerpart();
			nnparam->joiner->nearest_neighbor_aabb(p.first, p.second);
		}else{
			nnparam->joiner->nearest_neighbor(p.first, p.second);
		}
		if(p.second!=p.first){
			delete p.second;
		}
		delete p.first;
	}
	return NULL;
}

void SpatialJoin::nearest_neighbor_batch(vector<pair<Tile *, Tile *>> &tile_pairs, int num_threads, bool ispeed){
	struct nn_param param;
	for(pair<Tile *, Tile *> &p:tile_pairs){
		param.tile_queue.push(p);
	}
	param.joiner = this;
	param.ispeed = ispeed;
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
		if(p.second!=p.first){
			delete p.second;
		}
		delete p.first;
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


