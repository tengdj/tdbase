/*
 * WithinJoin.cpp
 *
 *  Created on: Sep 16, 2022
 *      Author: teng
 */

#include "SpatialJoin.h"

namespace hispeed{

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
	vector<pair<int, range>> candidate_ids;
	OctreeNode *tree = tile2->get_octree();
	size_t tile1_size = min(tile1->num_objects(), ctx.max_num_objects1);
	for(int i=0;i<tile1_size;i++){
		HiMesh_Wrapper *wrapper1 = tile1->get_mesh_wrapper(i);
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
					// must not within
					if(dist_vox.mindist>ctx.within_dist){
						continue;
					}
					// must be within
					if(dist_vox.maxdist<=ctx.within_dist){
						determined = true;
						break;
					}
					// the faces in those voxels need be further evaluated
					ci.voxel_pairs.push_back(voxel_pair(v1, v2, dist_vox));
					min_maxdist = min(min_maxdist, dist_vox.maxdist);
				}
				if(determined){
					break;
				}
			}

			// determined with the voxel evaluation
			if(determined){
				wrapper1->report_result(wrapper2);
				//delete ci;
				continue;
			}

			// otherwise, for further evaluation
			ci.distance = update_voxel_pair_list(ci.voxel_pairs, min_maxdist);
			// some voxel pairs need to be further evaluated
			if(ci.voxel_pairs.size()>0){
				ce->add_candidate(ci);
			}else{
				//delete ci;
			}
		}
		// save the candidate list
		if(ce->candidates.size()>0){
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
void SpatialJoin::within(query_context ctx){
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

		// now update the candidate list with the new latest information
		int index = 0;
		start = get_cur_time();
		for(auto ce_iter=candidates.begin();ce_iter!=candidates.end();){
			HiMesh_Wrapper *wrapper1 = (*ce_iter)->mesh_wrapper;
			//print_candidate_within(*ce_iter);
			for(auto ci_iter=(*ce_iter)->candidates.begin();ci_iter!=(*ce_iter)->candidates.end();){
				bool determined = false;
				HiMesh_Wrapper *wrapper2 = (ci_iter)->mesh_wrapper;
				if(ctx.use_aabb){
					range dist = (ci_iter)->distance;
					result_container res = ctx.results[index++];
					float hdist1 = wrapper1->getHausdorffDistance();
					float hdist2 = wrapper2->getHausdorffDistance();
					if(lod==ctx.highest_lod()){
						// now we have a precise distance
						dist.mindist = res.distance;
						dist.maxdist = res.distance;
					}else{
						dist.maxdist = std::min(dist.maxdist, res.distance);
						dist.mindist = std::max(dist.mindist, dist.maxdist-hdist1-hdist2);
						if(global_ctx.verbose>=1){
							log("%ld\t%ld:\t%.2f %.2f\t[%.2f, %.2f]->[%.2f, %.2f]",wrapper1->id, wrapper2->id,
									hdist1, hdist2,
									(ci_iter)->distance.mindist, (ci_iter)->distance.maxdist,
									dist.mindist, dist.maxdist);
						}
					}
					(ci_iter)->distance = dist;
					if(dist.maxdist<=ctx.within_dist){
						// the distance is close enough
						wrapper1->report_result(wrapper2);
						determined = true;
					}else if(dist.mindist > ctx.within_dist){
						// not possible
						determined = true;
					}
				}else{
					for(auto vp_iter = (ci_iter)->voxel_pairs.begin();vp_iter!=(ci_iter)->voxel_pairs.end();){
						result_container res = ctx.results[index++];
						//cout<<vp_iter->v1->num_triangles<<"  "<<res.p1<<" "<<vp_iter->v2->num_triangles<<" "<<res.p2<<" "<<res.distance<<endl;
						// update the distance
						if(!determined && vp_iter->v1->num_triangles>0&&vp_iter->v2->num_triangles>0){
							range dist = vp_iter->dist;
							float hdist1;
							float hdist2;
							float phdist1;
							float phdist2;
							if(global_ctx.hausdorf_level<2){
								hdist1 = wrapper1->getHausdorffDistance();
								hdist2 = wrapper2->getHausdorffDistance();
								phdist1 = wrapper1->getProxyHausdorffDistance();
								phdist2 = wrapper2->getProxyHausdorffDistance();
							}else{
								hdist1 = vp_iter->v1->getHausdorffDistance(res.p1);
								hdist2 = vp_iter->v2->getHausdorffDistance(res.p2);
								phdist1 = vp_iter->v1->getProxyHausdorffDistance(res.p1);
								phdist2 = vp_iter->v2->getProxyHausdorffDistance(res.p2);
							}

							if(lod==ctx.highest_lod()){
								// now we have a precise distance
								dist.mindist = res.distance;
								dist.maxdist = res.distance;
							}else{
								dist.maxdist = std::min(dist.maxdist, res.distance);
								if(global_ctx.hausdorf_level>0)
								{
									dist.mindist = std::max(dist.mindist, dist.maxdist-hdist1-hdist2);
									if(global_ctx.verbose>=1){
										log("%ld\t%ld:\t%.2f %.2f\t[%.2f, %.2f]->[%.2f, %.2f]",wrapper1->id, wrapper2->id,
												hdist1, hdist2,
												(ci_iter)->distance.mindist, (ci_iter)->distance.maxdist,
												dist.mindist, dist.maxdist);
									}
								}
							}
							vp_iter->dist = dist;
							// one voxel pair is close enough
							if(dist.maxdist<=ctx.within_dist){
								determined = true;
								wrapper1->report_result(wrapper2);
							}
						}
						// too far, should be removed from the voxel pair list
						if(vp_iter->dist.mindist>ctx.within_dist){
							(ci_iter)->voxel_pairs.erase(vp_iter);
						}else{
							vp_iter++;
						}
					}
				}

				if(determined || (ci_iter)->voxel_pairs.size()==0){
					// must closer than or farther than
					//delete *ci_iter;
					(*ce_iter)->candidates.erase(ci_iter);
				}else{
					ci_iter++;
				}
			}
			if((*ce_iter)->candidates.size()==0){
				delete (*ce_iter);
				candidates.erase(ce_iter);
			}else{
				//print_candidate_within(*ce_iter);
				ce_iter++;
			}
		}
		delete []ctx.results;
		ctx.updatelist_time += logt("updating the candidate lists",start);

		logt("evaluating with lod %d", iter_start, lod);
		log("");
	}
	ctx.overall_time = hispeed::get_time_elapsed(very_start, false);
	for(int i=0;i<ctx.tile1->num_objects();i++){
		ctx.result_count += ctx.tile1->get_mesh_wrapper(i)->results.size();
		for(int j=0;j<ctx.tile1->get_mesh_wrapper(i)->results.size();j++){
			//cout<<ctx.tile1->get_mesh_wrapper(i)->id<<"\t"<<ctx.tile1->get_mesh_wrapper(i)->results[j]->id<<endl;
		}
	}
	ctx.obj_count += min(ctx.tile1->num_objects(),global_ctx.max_num_objects1);
	global_ctx.merge(ctx);
}

}
