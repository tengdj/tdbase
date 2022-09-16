/*
 * WithinJoin.cpp
 *
 *  Created on: Sep 16, 2022
 *      Author: teng
 */

#include "SpatialJoin.h"

namespace hispeed{

inline void update_candidate_list_within(vector<candidate_entry> &candidates){
	for(vector<candidate_entry>::iterator it = candidates.begin();it!=candidates.end();){
		bool intersected = false;
		for(candidate_info &info:it->second){
			for(voxel_pair &vp:info.voxel_pairs){
				// if any voxel pair is ensured to be intersected
				if(vp.fulfill){
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

vector<candidate_entry> SpatialJoin::mbb_within(Tile *tile1, Tile *tile2, query_context &ctx){
	vector<candidate_entry> candidates;
	vector<pair<int, range>> candidate_ids;
	OctreeNode *tree = tile2->build_octree(20);
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

/*
 * the main function for getting the object within a specified distance
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
								vp.fulfill = true;
							}
						}
						index++;
					}
				}
			}
		}
		update_candidate_list_within(candidates);

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

}
