/*
 * WithinJoin.cpp
 *
 *  Created on: Sep 16, 2022
 *      Author: teng
 */

#include "SpatialJoin.h"

namespace tdbase{

inline void print_candidate_within(candidate_entry *cand){
	if(global_ctx.verbose>=1){
		printf("%ld (%ld candidates)\n", cand->mesh_wrapper->id, cand->candidates.size());
		for(int i=0;i<cand->candidates.size();i++){
			printf("%d:\t%ld\t%ld\n",i,cand->candidates[i].mesh_wrapper->id,cand->candidates[i].voxel_pairs.size());
			for(auto &vp:cand->candidates[i].voxel_pairs){
				printf("\t[%f,%f]\n",vp.dist.mindist,vp.dist.maxdist);
			}
		}
	}
}

vector<candidate_entry *> SpatialJoin::mbb_within(Tile *tile1, Tile *tile2, query_context &ctx){
	vector<candidate_entry *> candidates;
	OctreeNode *tree = tile2->get_octree();
	size_t tile1_size = min(tile1->num_objects(), ctx.max_num_objects1);
//#pragma omp parallel for
	for(int i=0;i<tile1_size;i++){
		vector<pair<int, range>> candidate_ids;
		HiMesh_Wrapper *wrapper1 = tile1->get_mesh_wrapper(i);
//		if(wrapper1->id!=7142){
//			continue;
//		}
		tree->query_within(&(wrapper1->box), candidate_ids, ctx.within_dist);
		if(candidate_ids.empty()){
			continue;
		}

		candidate_entry *ce = new candidate_entry(wrapper1);
		for(pair<int, range> &p:candidate_ids){
			HiMesh_Wrapper *wrapper2 = tile2->get_mesh_wrapper(p.first);
			candidate_info ci(wrapper2);
			bool determined = false;
			float min_maxdist = DBL_MAX;
			for(Voxel *v1:wrapper1->voxels){
				for(Voxel *v2:wrapper2->voxels){
					range dist_vox = v1->distance(*v2);
					// must not be within the specified distance
					if(dist_vox.mindist>ctx.within_dist){
						continue;
					}
					// must be within
					if(dist_vox.maxdist<=ctx.within_dist){
						determined = true;
						break;
					}
					// the facets in those voxels need be further evaluated
					ci.voxel_pairs.push_back(voxel_pair(v1, v2, dist_vox));
					min_maxdist = min(min_maxdist, dist_vox.maxdist);
				}
				if(determined){
					break;
				}
			}

			// determined with the voxel evaluation
			if(determined){
//#pragma omp critical
				ctx.report_result(wrapper1->id, wrapper2->id);
			}else if(ci.voxel_pairs.size()>0){
				// some voxel pairs need to be further evaluated
				ci.distance = update_voxel_pair_list(ci.voxel_pairs, min_maxdist);
				ce->add_candidate(ci);
			}
		}
		// save the candidate list
		if(ce->candidates.size()>0){
//#pragma omp critical
			candidates.push_back(ce);
		}else{
			delete ce;
		}
		candidate_ids.clear();
	}
	return candidates;
}

/*
 * the main function for getting the object within a specified distance
 *
 * */
void SpatialJoin::within_distance(query_context ctx){
	struct timeval start = get_cur_time();
	struct timeval very_start = get_cur_time();
	// filtering with MBBs to get the candidate list
	vector<candidate_entry *> candidates = mbb_within(ctx.tile1, ctx.tile2, ctx);
	ctx.index_time += get_time_elapsed(start, false);
	logt("comparing mbbs with %d candidate pairs", start, get_candidate_num(candidates));

	// now we start to get the distances with progressive level of details
	for(uint32_t lod:ctx.lods){
		ctx.cur_lod = lod;
		struct timeval iter_start = get_cur_time();
		const int pair_num = get_pair_num(candidates);
		if(pair_num==0){
			break;
		}
		size_t candidate_num = get_candidate_num(candidates);
		log("%ld polyhedron has %d candidates %d voxel pairs %.2f voxel pairs per candidate",
				candidates.size(), candidate_num, pair_num, (1.0*pair_num)/candidates.size());

		// do the computation
		calculate_distance(candidates, ctx);

		// now update the candidate list with the latest distance calculated with this LOD
		int index = 0;
		start = get_cur_time();
		for(auto ce_iter=candidates.begin();ce_iter!=candidates.end();){
			HiMesh_Wrapper *wrapper1 = (*ce_iter)->mesh_wrapper;
			//print_candidate_within(*ce_iter);
			for(auto ci_iter=(*ce_iter)->candidates.begin();ci_iter!=(*ce_iter)->candidates.end();){
				bool determined = false;
				HiMesh_Wrapper *wrapper2 = (ci_iter)->mesh_wrapper;

				if(ctx.use_aabb){
					// The AABB-based approach is disabled for progressive querying
					assert(lod==ctx.highest_lod());
					result_container res = ctx.tmp_results[index++];

					ci_iter->distance.mindist = res.distance;
					ci_iter->distance.maxdist = res.distance;

					if(res.distance<=ctx.within_dist){
						// the distance is close enough
						ctx.report_result(wrapper1->id, wrapper2->id);
						determined = true;
					}else if(res.distance > ctx.within_dist){
						// not possible
						determined = true;
					}
				}else{
					double vox_minmaxdist = DBL_MAX;
					for(auto &vp:ci_iter->voxel_pairs){
						result_container res = ctx.tmp_results[index++];

						// skip evaluating this pair in this round
						// no geometric computation is conducted in this round
						if(vp.has_empty_voxel()){
							if(global_ctx.verbose>=1) {
								log("%ld(%d)\t%ld(%d):"
										"[%.2f, %.2f]"
										" has empty voxel",
										wrapper1->id, vp.v1->id,
										wrapper2->id, vp.v2->id,
										vp.dist.mindist, vp.dist.maxdist);
							}
							continue;
						}

						// update the distance
						range dist = vp.dist;
						if(lod==ctx.highest_lod()){
							// now we have a precise distance
							dist.mindist = res.distance;
							dist.maxdist = res.distance;
						} else {
							dist.mindist = std::max(dist.mindist, res.min_dist);
							dist.maxdist = std::min(dist.maxdist, res.max_dist);
						}
						//dist.maxdist = std::min(dist.maxdist, res.distance);

						if(global_ctx.verbose>=1) {
							log("%ld(%d)\t%ld(%d):"
									"[%.2f, %.2f]"
									"->[%.2f, %.2f]",
									wrapper1->id, vp.v1->id,
									wrapper2->id, vp.v2->id,
									vp.dist.mindist, vp.dist.maxdist,
									dist.mindist, dist.maxdist);
						}
						vp.dist = dist;
						vox_minmaxdist = min(vox_minmaxdist, (double)vp.dist.maxdist);
					}// end for voxel pair iteration

					assert(ci_iter->voxel_pairs.size()>0);
					// update the object-level distance with the voxel pairs
					ci_iter->distance = update_voxel_pair_list(ci_iter->voxel_pairs, vox_minmaxdist,lod!=ctx.highest_lod());
				}// end AABB or progressive option

				if(ci_iter->distance.mindist > ctx.within_dist){
					// the object must farther than the target
					(*ce_iter)->candidates.erase(ci_iter);
				}else if(ci_iter->distance.maxdist <= ctx.within_dist){
					// the object must closer than the target
					ctx.report_result(wrapper1->id, wrapper2->id);
					(*ce_iter)->candidates.erase(ci_iter);
				}else{
					ci_iter++;
				}
			}
			// all candidates for one object are eliminated, remove such entry
			if((*ce_iter)->candidates.size()==0){
				delete (*ce_iter);
				candidates.erase(ce_iter);
			}else{
				//print_candidate_within(*ce_iter);
				ce_iter++;
			}
		}
		delete []ctx.tmp_results;
		ctx.updatelist_time += logt("updating the candidate lists",start);

		logt("evaluating with lod %d", iter_start, lod);
		log("");
	}
	ctx.overall_time = tdbase::get_time_elapsed(very_start, false);
	ctx.obj_count += min(ctx.tile1->num_objects(),global_ctx.max_num_objects1);
	global_ctx.merge(ctx);
}

}
