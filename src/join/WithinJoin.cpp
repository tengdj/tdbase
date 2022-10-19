/*
 * WithinJoin.cpp
 *
 *  Created on: Sep 16, 2022
 *      Author: teng
 */

#include "SpatialJoin.h"

namespace hispeed{

inline void print_candidate_within(candidate_entry &cand){
	if(global_ctx.verbose){
		printf("%ld (%ld candidates)\n", cand.mesh_wrapper->id, cand.candidates.size());
		for(int i=0;i<cand.candidates.size();i++){
			printf("%d:\t%ld\t%ld\n",i,cand.candidates[i].mesh_wrapper->id,cand.candidates[i].voxel_pairs.size());
			for(auto &vp:cand.candidates[i].voxel_pairs){
				printf("\t[%f,%f]\n",vp.dist.mindist,vp.dist.maxdist);
			}
		}
	}
}

vector<candidate_entry> SpatialJoin::mbb_within(Tile *tile1, Tile *tile2, query_context &ctx){
	vector<candidate_entry> candidates;
	vector<pair<int, range>> candidate_ids;
	OctreeNode *tree = tile2->build_octree(20);
	size_t tile1_size = min(tile1->num_objects(), ctx.max_num_objects1);
	for(int i=0;i<tile1_size;i++){
		vector<candidate_info> candidate_list;
		HiMesh_Wrapper *wrapper1 = tile1->get_mesh_wrapper(i);
		tree->query_within(&(wrapper1->box), candidate_ids, ctx.max_dist);
		if(candidate_ids.empty()){
			continue;
		}
		for(pair<int, range> &p:candidate_ids){
			HiMesh_Wrapper *wrapper2 = tile2->get_mesh_wrapper(p.first);
			candidate_info ci;
			bool determined = false;
			float min_maxdist = DBL_MAX;
			for(Voxel *v1:wrapper1->voxels){
				for(Voxel *v2:wrapper2->voxels){
					range dist_vox = v1->distance(*v2);
					// must not within
					if(dist_vox.mindist>ctx.max_dist){
						continue;
					}
					// must be within
					if(dist_vox.maxdist<=ctx.max_dist){
						determined = true;
						break;
					}
					// the faces in those voxels need be further evaluated
					ci.voxel_pairs.push_back(voxel_pair(v1, v2, dist_vox));
					min_maxdist = min(min_maxdist, dist_vox.maxdist);
				}
				if(determined){
					wrapper1->report_result(wrapper2);
					ci.voxel_pairs.clear();
					break;
				}
			}
			ci.distance = update_voxel_pair_list(ci.voxel_pairs, min_maxdist);
			// some voxel pairs need to be further evaluated
			if(ci.voxel_pairs.size()>0){
				ci.mesh_wrapper = wrapper2;
				candidate_list.push_back(ci);
			}
		}
		// save the candidate list
		if(candidate_list.size()>0){
			candidates.push_back(candidate_entry(wrapper1, candidate_list));
		}
		candidate_ids.clear();
	}
	delete tree;
	return candidates;
}

/*
 * the main function for getting the object within a specified distance
 *
 * */
void SpatialJoin::within(query_context ctx){
	struct timeval start = get_cur_time();
	struct timeval very_start = get_cur_time();

	// filtering with MBBs to get the candidate list
	vector<candidate_entry> candidates = mbb_within(ctx.tile1, ctx.tile2, ctx);
	ctx.index_time += get_time_elapsed(start, false);
	logt("comparing mbbs", start);

	// now we start to get the distances with progressive level of details
	for(int lod:ctx.lods){
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

		// now update the candidate list with the new latest information
		int index = 0;
		start = get_cur_time();
		for(auto ce_iter=candidates.begin();ce_iter!=candidates.end();){
			HiMesh_Wrapper *wrapper1 = ce_iter->mesh_wrapper;
			//print_candidate_within(*ce_iter);
			for(auto ci_iter=ce_iter->candidates.begin();ci_iter!=ce_iter->candidates.end();){
				bool determined = false;
				HiMesh_Wrapper *wrapper2 = ci_iter->mesh_wrapper;
				if(ctx.use_aabb){
					range dist = ci_iter->distance;
					if(lod==ctx.highest_lod()){
						// now we have a precise distance
						dist.mindist = ctx.distance[index];
						dist.maxdist = ctx.distance[index];
					}else{
						dist.maxdist = std::min(dist.maxdist, ctx.distance[index]);
						dist.mindist = std::max(dist.mindist, dist.maxdist-wrapper1->mesh->getmaximumCut()-wrapper2->mesh->getmaximumCut());
						if(global_ctx.verbose){
							log("%ld\t%ld:\t%.2f %.2f\t[%.2f, %.2f]->[%.2f, %.2f]",wrapper1->id, wrapper2->id,
									wrapper1->mesh->getmaximumCut(), wrapper2->mesh->getmaximumCut(),
									ci_iter->distance.mindist, ci_iter->distance.maxdist,
									dist.mindist, dist.maxdist);
						}
					}
					ci_iter->distance = dist;
					if(dist.maxdist<=ctx.max_dist){
						// the distance is close enough
						wrapper1->report_result(wrapper2);
						determined = true;
					}else if(dist.mindist > ctx.max_dist){
						// not possible
						determined = true;
					}
					index++;
				}else{
					for(auto vp_iter = ci_iter->voxel_pairs.begin();vp_iter!=ci_iter->voxel_pairs.end();){
						// update the distance
						if(!determined && vp_iter->v1->size[lod]>0&&vp_iter->v2->size[lod]>0){
							range dist = vp_iter->dist;
							if(lod==ctx.highest_lod()){
								// now we have a precise distance
								dist.mindist = ctx.distance[index];
								dist.maxdist = ctx.distance[index];
							}else{
								dist.maxdist = std::min(dist.maxdist, ctx.distance[index]);
								dist.mindist = std::max(dist.mindist, dist.maxdist-wrapper1->mesh->getmaximumCut()-wrapper2->mesh->getmaximumCut());
								if(global_ctx.verbose){
									log("%ld\t%ld:\t%.2f %.2f\t[%.2f, %.2f]->[%.2f, %.2f]",wrapper1->id, wrapper2->id,
											wrapper1->mesh->getmaximumCut(), wrapper2->mesh->getmaximumCut(),
											ci_iter->distance.mindist, ci_iter->distance.maxdist,
											dist.mindist, dist.maxdist);
								}
							}
							vp_iter->dist = dist;
							// one voxel pair is close enough
							if(dist.maxdist<=ctx.max_dist){
								determined = true;
								wrapper1->report_result(wrapper2);
							}
						}
						// too far, should be removed from the voxel pair list
						if(vp_iter->dist.mindist>ctx.max_dist){
							ci_iter->voxel_pairs.erase(vp_iter);
						}else{
							vp_iter++;
						}
						index++;
					}
				}
				if(determined || ci_iter->voxel_pairs.size()==0){
					// must closer than or farther than
					ce_iter->candidates.erase(ci_iter);
				}else{
					ci_iter++;
				}
			}
			if(ce_iter->candidates.size()==0){
				candidates.erase(ce_iter);
			}else{
				//print_candidate_within(*ce_iter);
				ce_iter++;
			}
		}
		delete []ctx.distance;
		ctx.updatelist_time += logt("updating the candidate lists",start);

		logt("evaluating with lod %d", iter_start, lod);
		log("");
	}
	ctx.overall_time = hispeed::get_time_elapsed(very_start, false);
	for(int i=0;i<ctx.tile1->num_objects();i++){
		ctx.result_count += ctx.tile1->get_mesh_wrapper(i)->results.size();
	}
	ctx.obj_count += min(ctx.tile1->num_objects(),global_ctx.max_num_objects1);
	global_ctx.merge(ctx);
}

}
