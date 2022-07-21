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


/*facilitate functions*/


inline bool update_voxel_pair_list(vector<voxel_pair> &voxel_pairs, range &d){
	int voxel_pair_size = voxel_pairs.size();
	for(int j=0;j<voxel_pair_size;){
		if(d>voxel_pairs[j].dist){
			return false;
		}else if(d<voxel_pairs[j].dist){
			// evict this dequalified voxel pairs
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
			// one candidate in the list cannot be the mindist
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


/*
 *
 * functions for checking indexes
 *
 */

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
					if(v1->intersect(*v2)){
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



vector<candidate_entry> SpatialJoin::mbb_nn(Tile *tile1, Tile *tile2, query_context &ctx){
	vector<candidate_entry> candidates;
	vector<pair<int, range>> candidate_ids;
	OctreeNode *tree = tile2->build_octree(400);
	for(int i=0;i<tile1->num_objects();i++){
		// for each object
		//1. use the distance between the mbbs of objects as a
		//	 filter to retrieve candidate objects
		vector<candidate_info> candidate_list;
		HiMesh_Wrapper *wrapper1 = tile1->get_mesh_wrapper(i);
		float min_maxdistance = DBL_MAX;
		tree->query_nn(&(wrapper1->box), candidate_ids, min_maxdistance);
		if(candidate_ids.empty()){
			continue;
		}
		//cout<<i<<" tengdj "<<candidate_ids.size()<<endl;
		std::sort(candidate_ids.begin(), candidate_ids.end(), compare_pair);
		int former = -1;
		for(pair<int, range> &p:candidate_ids){
			if(p.first==former){
				// duplicate
				continue;
			}
			former = p.first;

			HiMesh_Wrapper *wrapper2 = tile2->get_mesh_wrapper(p.first);
			//2. we further go through the voxels in two objects to shrink
			// 	 the candidate list in a fine grained

			//2.1 first go over the candidate list with current distance range
			//
			if(update_candidate_list(candidate_list, p.second)){
				candidate_info ci;
				ci.mesh_wrapper = wrapper2;
				for(Voxel *v1:wrapper1->voxels){
					for(Voxel *v2:wrapper2->voxels){
						range tmpd = v1->distance(*v2);
						// no voxel pair in the list is nearer
						if(update_voxel_pair_list(ci.voxel_pairs, tmpd) &&
						   update_candidate_list(candidate_list, tmpd)){
							ci.voxel_pairs.push_back(voxel_pair(v1, v2, tmpd));
						}
					}
				}
				// some voxel pairs need be further evaluated
				if(ci.voxel_pairs.size()>0){
					ci.distance = ci.voxel_pairs[0].dist;
					for(voxel_pair &p:ci.voxel_pairs){
						ci.distance.update(p.dist);
					}
					assert(ci.distance.mindist<=ci.distance.maxdist);
					candidate_list.push_back(ci);
				}
			}

		}
		// save the candidate list
		if(candidate_list.size()>0){
			candidates.push_back(candidate_entry(wrapper1, candidate_list));
		}
		candidate_ids.clear();
	}
	return candidates;
}

vector<candidate_entry> SpatialJoin::mbb_within(Tile *tile1, Tile *tile2, query_context &ctx){
	vector<candidate_entry> candidates;
	vector<pair<int, range>> candidate_ids;
	OctreeNode *tree = tile2->build_octree(400);
	for(int i=0;i<tile1->num_objects();i++){
		vector<candidate_info> candidate_list;
		HiMesh_Wrapper *wrapper1 = tile1->get_mesh_wrapper(i);
		tree->query_within(&(wrapper1->box), candidate_ids, ctx.max_dist);
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

			candidate_info ci;
			bool determined = false;
			for(Voxel *v1:wrapper1->voxels){
				for(Voxel *v2:wrapper2->voxels){
					range tmpd = v1->distance(*v2);
					// must not within
					if(tmpd.mindist>ctx.max_dist){
						continue;
					}
					// must be within
					if(tmpd.maxdist<=ctx.max_dist){
						determined = true;
						//todo report here
					}
					if(determined){
						break;
					}
					// no voxel pair in the list is nearer
					ci.voxel_pairs.push_back(voxel_pair(v1, v2, tmpd));
				}
				if(determined){
					ci.voxel_pairs.clear();
					break;
				}
			}
			// some voxel pairs need be further evaluated
			if(ci.voxel_pairs.size()>0){
				ci.mesh_wrapper = wrapper2;
				candidate_list.push_back(ci);
			}

			former = p.first;
		}
		// save the candidate list
		if(candidate_list.size()>0){
			candidates.push_back(candidate_entry(wrapper1, candidate_list));
		}
		candidate_ids.clear();
	}
	return candidates;
}

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
		for(candidate_entry &c:candidates){
			// the nearest neighbor is found
			if(ctx.query_type=="nn"&&c.second.size()<=1){
				continue;
			}
			for(candidate_info &info:c.second){
				info.mesh_wrapper->mesh->get_aabb_tree_segment();
			}
		}
		ctx.packing_time += hispeed::get_time_elapsed(start, false);
		logt("building aabb tree", start);

		int index = 0;
		for(candidate_entry c:candidates){
			// the nearest neighbor is found
			if(c.second.size()<=1){
				continue;
			}
			vector<Point> vertices;
			c.first->mesh->get_vertices(vertices);
			for(candidate_info info:c.second){
				double min_dist = DBL_MAX;
				HiMesh_Wrapper *wrapper2 = info.mesh_wrapper;
				for(Point &p:vertices){
					FT sqd = wrapper2->mesh->get_aabb_tree_segment()->squared_distance(p);
					double dist = (double)CGAL::to_double(sqd);
//					Point tp = wrapper2->mesh->get_aabb_tree()->closest_point(p);
//					double  dist = distance(p,tp);
					if(min_dist>dist){
						min_dist = dist;
					}
				}
				assert(info.voxel_pairs.size()==1);
				distances[index++] = min_dist;
				wrapper2->mesh->clear_aabb_tree();
			}// end for distance_candiate list
			vertices.clear();
		}// end for candidates

	}else{

		map<Voxel *, std::pair<uint, uint>> voxel_map;
		uint segment_num = 0;

		for(candidate_entry &c:candidates){
			HiMesh_Wrapper *wrapper1 = c.first;
			for(candidate_info &info:c.second){
				for(voxel_pair &vp:info.voxel_pairs){
					assert(vp.v1&&vp.v2);
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
				}
			}
		}

		assert(segment_num>0 && "there should be segments in voxels need be calculated");

		// now we allocate the space and store the data in a buffer
		float *data = new float[6*segment_num];
		uint *offset_size = new uint[4*pair_num];

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
		ctx.packing_time += hispeed::get_time_elapsed(start, false);
		logt("organizing data", start);
		geometry_param gp;
		gp.data = data;
		gp.pair_num = pair_num;
		gp.offset_size = offset_size;
		gp.distances = distances;
		gp.data_size = segment_num;
		computer->get_distance(gp);
		delete data;
		delete offset_size;

	}

	return distances;
}


/*
 * the main function for getting the objecst within a specified distance
 *
 * */
void SpatialJoin::within(Tile *tile1, Tile *tile2, query_context ctx){
	struct timeval start = get_cur_time();
	struct timeval very_start = get_cur_time();

	// filtering with MBBs to get the candidate list
	vector<candidate_entry> candidates = mbb_within(tile1, tile2, ctx);
	ctx.index_time += get_time_elapsed(start, false);
	logt("comparing mbbs", start);

	// now we start to get the distances with progressive level of details
	for(int lod:ctx.lods){
		struct timeval iter_start = get_cur_time();
		const int pair_num = get_pair_num(candidates);
		if(pair_num==0){
			break;
		}
		size_t candidate_num = get_candidate_num(candidates);
		log("%ld polyhedron has %d candidates %f voxel pairs per candidate", candidates.size(), candidate_num, (1.0*pair_num)/candidates.size());
		// retrieve the necessary meshes
		size_t segment_pair_num = 0;


		int o1_counter = 0;
		int o2_counter = 0;
		int o1_counter_nocache = 0;
		int o2_counter_nocache = 0;

		for(candidate_entry &c:candidates){
			HiMesh_Wrapper *wrapper1 = c.first;
			o1_counter_nocache++;
			for(candidate_info &info:c.second){
				HiMesh_Wrapper *wrapper2 = info.mesh_wrapper;
				for(voxel_pair &vp:info.voxel_pairs){
					assert(vp.v1&&vp.v2);
					// not filled yet
					if(vp.v1->data.find(lod)==vp.v1->data.end()){
						// ensure the mesh is extracted
						tile1->decode_to(wrapper1->id, lod);
						wrapper1->fill_voxels(DT_Segment);
						o1_counter++;
					}
					if(vp.v2->data.find(lod)==vp.v2->data.end()){
						tile2->decode_to(wrapper2->id, lod);
						wrapper2->fill_voxels(DT_Segment);
						o2_counter++;
					}
					o2_counter_nocache++;

					segment_pair_num += vp.v1->size[lod]*vp.v2->size[lod];
				}// end for voxel_pairs
			}// end for distance_candiate list
		}// end for candidates
		ctx.decode_time += hispeed::get_time_elapsed(start, false);
		logt("decoded %ld segment pairs for lod %d", start, segment_pair_num, lod);
		if(segment_pair_num==0){
			log("no segments is filled in this round");
			continue;
		}
		tile1->reset_time();
		log("%d,%d,%d,%d,%d",lod,o1_counter,o1_counter_nocache,o2_counter,o2_counter_nocache);


		float *distances = this->calculate_distance(candidates, ctx, lod);

		ctx.computation_time += hispeed::get_time_elapsed(start, false);
		logt("get distance", start);

		// now update the candidate list with the new distance information
		int index = 0;
		for(candidate_entry &ce:candidates){
			bool determined = false;
			for(candidate_info &c:ce.second){
				if(!determined){
					for(voxel_pair &vp:c.voxel_pairs){
						// update the distance
						if(!determined&&vp.v1->size[lod]>0&&vp.v2->size[lod]>0){
							range dist = vp.dist;

							if(lod==ctx.highest_lod()){
								// now we have a precise distance
								dist.mindist = distances[index];
								dist.maxdist = distances[index];
							}else{
								dist.maxdist = std::min(dist.maxdist, distances[index]);
							}
							vp.dist = dist;
							if(dist.maxdist<ctx.max_dist){
								determined = true;
								vp.intersect = true;
							}
						}
						index++;
					}
				}
			}
		}
		hispeed::update_candidate_list_intersect(candidates);

		ctx.updatelist_time += hispeed::get_time_elapsed(start, false);
		logt("update candidate list", start);

		delete distances;
		logt("current iteration", iter_start);

	}
	ctx.overall_time = hispeed::get_time_elapsed(very_start, false);
	pthread_mutex_lock(&g_lock);
	global_ctx.merge(ctx);
	pthread_mutex_unlock(&g_lock);

}




/*
 * the main function for getting the nearest neighbor
 *
 * */

int cc = 0;
void SpatialJoin::nearest_neighbor(Tile *tile1, Tile *tile2, query_context ctx){
	struct timeval start = get_cur_time();
	struct timeval very_start = get_cur_time();

	// filtering with MBBs to get the candidate list
	vector<candidate_entry> candidates = mbb_nn(tile1, tile2, ctx);
	ctx.index_time += get_time_elapsed(start, false);
	for(vector<candidate_entry>::iterator it=candidates.begin();it!=candidates.end();){
		if(it->second.size()<=1){
			it = candidates.erase(it);
		}else{
			it++;
		}
	}
	logt("comparing mbbs", start);

	// now we start to get the distances with progressive level of details
	for(int lod:ctx.lods){
		struct timeval iter_start = get_cur_time();
		const int pair_num = get_pair_num(candidates);
		if(pair_num==0){
			break;
		}
		size_t candidate_num = get_candidate_num(candidates);
		log("%ld polyhedron has %d candidates %d voxel pairs %f voxel pairs per candidate",
				candidates.size(), candidate_num, pair_num, (1.0*pair_num)/candidates.size());
		// retrieve the necessary meshes
		size_t segment_pair_num = 0;

		int o1_counter = 0;
		int o2_counter = 0;
		int o2_counter_nocache = 0;
		int o1_counter_nocache = 0;

		for(candidate_entry &c:candidates){
			// the nearest neighbor is found
			if(c.second.size()<=1){
				continue;
			}
			o1_counter_nocache++;
			HiMesh_Wrapper *wrapper1 = c.first;
			for(candidate_info &info:c.second){
				for(voxel_pair &vp:info.voxel_pairs){
					assert(vp.v1&&vp.v2);
					// not filled yet
					if(vp.v1->data.find(lod)==vp.v1->data.end()){
						// ensure the mesh is extracted
						tile1->decode_to(wrapper1->id, lod);
						wrapper1->fill_voxels(DT_Segment);
						o1_counter++;
					}

					if(vp.v2->data.find(lod)==vp.v2->data.end()){
						tile2->decode_to(info.mesh_wrapper->id, lod);
						info.mesh_wrapper->fill_voxels(DT_Segment);
						o2_counter++;
					}
					o2_counter_nocache++;

					segment_pair_num += vp.v1->size[lod]*vp.v2->size[lod];
				}// end for voxel_pairs
			}// end for distance_candiate list
		}// end for candidates
		ctx.decode_time += hispeed::get_time_elapsed(start, false);
		logt("decoded with %ld segment pairs for lod %d", start, segment_pair_num, lod);
		if(segment_pair_num==0){
			log("no segments is filled in this round");
			continue;
		}

		tile1->reset_time();

		log("%d,%d,%d,%d,%d",lod,o1_counter,o1_counter_nocache,o2_counter,o2_counter_nocache);

		float *distances = this->calculate_distance(candidates, ctx, lod);

		ctx.computation_time += hispeed::get_time_elapsed(start, false);
		logt("get distance", start);

		// now update the distance range with the new distances
		int index = 0;
		for(candidate_entry &c:candidates){
			range min_candidate;
			min_candidate.maxdist = DBL_MAX;
			for(candidate_info &ci:c.second){
				for(voxel_pair &vp:ci.voxel_pairs){
					// update the distance
					if(vp.v1->size[lod]>0&&vp.v2->size[lod]>0){
						range dist = vp.dist;
						if(lod==ctx.highest_lod()){
							// now we have a precise distance
							dist.mindist = distances[index];
							dist.maxdist = distances[index];
						}else{
							dist.maxdist = std::min(dist.maxdist, distances[index]);
						}
						vp.dist = dist;
						if(min_candidate.maxdist>dist.maxdist){
							min_candidate = dist;
						}
					}
					index++;
				}
			}
			if(min_candidate.maxdist>ctx.max_nearest_distance){
				ctx.max_nearest_distance = min_candidate.maxdist;
			}
			update_candidate_list(c.second, min_candidate);
		}
		for(int i=0;i<candidates.size();i++){
			int target = 0;
			double min_dist = DBL_MAX;
			for(int j=0;j<candidates[i].second.size();j++){
				if(candidates[i].second[j].distance.maxdist<min_dist){
					min_dist = candidates[i].second[j].distance.maxdist;
					target = candidates[i].second[j].mesh_wrapper->id;
				}
			}
			//printf("%d\t%d\n",i,target);
		}
		for(vector<candidate_entry>::iterator it=candidates.begin();it!=candidates.end();){
			if(it->second.size()<=1){
				it = candidates.erase(it);
			}else{
				it++;
			}
		}
		ctx.updatelist_time += hispeed::get_time_elapsed(start, false);
		logt("update candidate list", start);

		delete distances;
		logt("current iteration", iter_start);
	}
	ctx.overall_time = hispeed::get_time_elapsed(very_start, false);
	pthread_mutex_lock(&g_lock);
	global_ctx.merge(ctx);
	pthread_mutex_unlock(&g_lock);
}



/*
 *
 * for doing intersection
 *
 * */




/*
 * the main function for detecting the intersection
 * relationship among polyhedra in the tile
 *
 * */
void SpatialJoin::intersect(Tile *tile1, Tile *tile2, query_context ctx){
	struct timeval start = get_cur_time();
	struct timeval very_start = get_cur_time();

	// filtering with MBBs to get the candidate list
	vector<candidate_entry> candidates = mbb_intersect(tile1, tile2);
	ctx.index_time += hispeed::get_time_elapsed(start,false);
	logt("comparing mbbs", start);

	// now we start to ensure the intersection with progressive level of details
	size_t triangle_pair_num = 0;

	for(int lod:ctx.lods){
		struct timeval iter_start = start;
		size_t pair_num = get_pair_num(candidates);
		if(pair_num==0){
			break;
		}
		log("%ld polyhedron has %ld candidates", candidates.size(), pair_num);
		// retrieve the necessary meshes
		map<Voxel *, std::pair<uint, uint>> voxel_map;
		size_t triangle_num = 0;

		int o1_counter = 0;
		int o2_counter = 0;
		int o2_counter_nocache = 0;
		int o1_counter_nocache = 0;
		for(candidate_entry c:candidates){
			HiMesh_Wrapper *wrapper1 = c.first;
			o1_counter_nocache++;
			for(candidate_info info:c.second){
				HiMesh_Wrapper *wrapper2 = info.mesh_wrapper;
				for(voxel_pair vp:info.voxel_pairs){
					// not filled yet
					if(vp.v1->data.find(lod)==vp.v1->data.end())
					{
						tile1->decode_to(wrapper1->id, lod);
						wrapper1->fill_voxels(DT_Triangle);
						o1_counter++;
					}
					if(vp.v2->data.find(lod)==vp.v2->data.end())
					{
						tile2->decode_to(wrapper2->id, lod);
						wrapper2->fill_voxels(DT_Triangle);
						o2_counter++;
					}
					o2_counter_nocache++;

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
		ctx.decode_time += hispeed::get_time_elapsed(start, false);
		logt("decoded %ld voxels with %ld triangles %ld pairs for lod %d",
				start, voxel_map.size(), triangle_num, triangle_pair_num, lod);

		//log("%d,%d,%d,%d,%d",lod,o1_counter,o1_counter_nocache,o2_counter,o2_counter_nocache);

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
		uint *intersect_status = new uint[pair_num];
		for(int i=0;i<pair_num;i++){
			intersect_status[i] = 0;
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
		ctx.packing_time += hispeed::get_time_elapsed(start, false);
		logt("organizing data", start);
		geometry_param gp;
		gp.data = data;
		gp.pair_num = pair_num;
		gp.offset_size = offset_size;
		gp.intersect = intersect_status;
		gp.data_size = triangle_num;
		computer->get_intersect(gp);
		ctx.computation_time += hispeed::get_time_elapsed(start, false);
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
		ctx.updatelist_time += hispeed::get_time_elapsed(start, false);
		logt("update candidate list", start);

		delete []data;
		delete []offset_size;
		delete []intersect_status;
		voxel_map.clear();

		logt("current iteration", iter_start);
	}
	ctx.overall_time = hispeed::get_time_elapsed(very_start, false);

	pthread_mutex_lock(&g_lock);
	global_ctx.merge(ctx);
	pthread_mutex_unlock(&g_lock);

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

void SpatialJoin::join(vector<pair<Tile *, Tile *>> &tile_pairs, query_context &ctx){
	struct nn_param param;
	for(pair<Tile *, Tile *> &p:tile_pairs){
		param.tile_queue.push(p);
	}
	param.joiner = this;
	param.ctx = ctx;
	pthread_t threads[ctx.num_thread];
	for(int i=0;i<ctx.num_thread;i++){
		pthread_create(&threads[i], NULL, join_unit, (void *)&param);
	}
	while(!param.tile_queue.empty()){
		usleep(10);
	}
	for(int i = 0; i < ctx.num_thread; i++){
		void *status;
		pthread_join(threads[i], &status);
	}

}


}


