/*
 * WithinJoin.cpp
 *
 *  Created on: Sep 16, 2022
 *      Author: teng
 */

#include "SpatialJoin.h"

namespace tdbase{

void DWithinJoin::index_retrieval(Tile *tile1, Tile *tile2, query_context &ctx){
	struct timeval start = get_cur_time();

	OctreeNode *tree = tile2->get_octree();

#pragma omp parallel for
	for(int i=0;i<tile1->num_objects();i++){
		vector<pair<int, range>> candidate_ids;
		HiMesh_Wrapper *wrapper1 = tile1->get_mesh_wrapper(i);

		if(config.specify_object!=-1&&config.specify_object!=wrapper1->id){//for single query
			continue;
		}

		tree->query_within(&(wrapper1->box), candidate_ids, config.within_dist);
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
					if(dist_vox.mindist>config.within_dist){
						continue;
					}
					// must be within
					if(dist_vox.maxdist<=config.within_dist){
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
#pragma omp critical
				ctx.report_result(wrapper1->id, wrapper2->id);
			}else if(ci.voxel_pairs.size()>0){
				// some voxel pairs need to be further evaluated
				ci.distance = update_voxel_pair_list(ci.voxel_pairs, min_maxdist);
				ce->add_candidate(ci);
			}
		}
		// save the candidate list
		if(ce->candidates.size()>0){
#pragma omp critical
			ctx.candidates.push_back(ce);
		}else{
			delete ce;
		}
		candidate_ids.clear();
	}

	ctx.index_time += logt("index retrieving", start);
	//evaluate_candidate_lists(ctx);
}

void DWithinJoin::evaluate_candidate_lists(query_context &ctx){
	struct timeval start = get_cur_time();
	update_distance_ranges(ctx);

	// evaluate the candidate lists for all the queried objects
	for(auto ce_iter=ctx.candidates.begin();ce_iter!=ctx.candidates.end();){
		HiMesh_Wrapper *wrapper1 = (*ce_iter)->mesh_wrapper;

		for(auto ci_iter=(*ce_iter)->candidates.begin();ci_iter!=(*ce_iter)->candidates.end();){
			HiMesh_Wrapper *wrapper2 = (ci_iter)->mesh_wrapper;
			if(ci_iter->distance.mindist > config.within_dist){
				// the object must farther than the target
				(*ce_iter)->candidates.erase(ci_iter);
			}else if(ci_iter->distance.maxdist <= config.within_dist){
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
			ctx.candidates.erase(ce_iter);
		}else{
			//print_candidate_within(*ce_iter);
			ce_iter++;
		}
	}

	ctx.updatelist_time += logt("updating the candidate lists",start);
}

}
