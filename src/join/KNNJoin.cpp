/*
 * KNNJoin.cpp
 *
 *  Created on: Sep 16, 2022
 *      Author: teng
 */

#include "SpatialJoin.h"

namespace tdbase{

void KNNJoin::index_retrieval(Tile *tile1, Tile *tile2, query_context &ctx){
	struct timeval start = get_cur_time();

	OctreeNode *tree = tile2->get_octree();

#pragma omp parallel for
	for(int i=0;i<tile1->num_objects();i++){
		vector<pair<int, range>> candidate_ids;
		// for each object
		//1. use the distance between the mbbs of objects as a
		//	 filter to retrieve candidate objects
		HiMesh_Wrapper *wrapper1 = tile1->get_mesh_wrapper(i);

		if(config.specify_object!=-1&&config.specify_object!=wrapper1->id){// for single query
			continue;
		}

		float min_maxdistance = DBL_MAX;
		tree->query_knn(&(wrapper1->box), candidate_ids, min_maxdistance, config.knn);
		assert(candidate_ids.size()>=config.knn);

		// result determined with only evaluating the MBBs
		if(candidate_ids.size() == config.knn){
			for(pair<int, range> &p:candidate_ids){
#pragma omp critical
				ctx.report_result(wrapper1->id, tile2->get_mesh_wrapper(p.first)->id);
			}
			candidate_ids.clear();
			continue;
		}

		candidate_entry *ce = new candidate_entry(wrapper1);

		//2. we further go through the voxels in two objects to shrink
		// 	 the candidate list in a finer granularity
		for(pair<int, range> &p:candidate_ids){
			HiMesh_Wrapper *wrapper2 = tile2->get_mesh_wrapper(p.first);

			candidate_info ci(wrapper2);
			float min_maxdist = DBL_MAX;

			for(Voxel *v1:wrapper1->voxels){
				for(Voxel *v2:wrapper2->voxels){
					range dist_vox = v1->distance(*v2);
					if(dist_vox.mindist>=min_maxdist){
						continue;
					}
					// wait for later evaluation
					ci.voxel_pairs.push_back(voxel_pair(v1, v2, dist_vox));
					min_maxdist = min(min_maxdist, dist_vox.maxdist);
				}
			}

			// form the distance range of objects with the evaluations of voxel pairs
			// should be updated after evaluating the voxels, as some pair maybe disqualified
			ci.distance = update_voxel_pair_list(ci.voxel_pairs, min_maxdist);
			assert(ci.voxel_pairs.size()>0);
			assert(ci.distance.valid());
			ce->add_candidate(ci);
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
	// the candidates list need be evaluated after checking with the mbb
	// some queries might be answered with only querying the index
	evaluate_candidate_lists(ctx);
}

/*
 *
 * the function for evaluating each queried object and its corresponding candidates
 *
 * */
void KNNJoin::evaluate_candidate_lists(query_context &ctx){
	struct timeval start = get_cur_time();

	update_distance_ranges(ctx);

#pragma omp parallel for
	for (candidate_entry *cand:ctx.candidates) {
		HiMesh_Wrapper *target = cand->mesh_wrapper;
		int list_size = cand->candidates.size();

		for(int i=0;i<list_size && config.knn>cand->candidate_confirmed;){
			int sure_closer = 0;
			int maybe_closer = 0;

			for(int j=0;j<cand->candidates.size();j++){
				if(i==j){
					continue;
				}
				// count how many candidates that are surely closer than this one
				if(cand->candidates[i].distance>=cand->candidates[j].distance) {
					sure_closer++;
				}
				// count how many candidates that are possibly closer than this one
				if(!(cand->candidates[i].distance<=cand->candidates[j].distance)) {
					maybe_closer++;
				}
			}
			int cand_left = config.knn-cand->candidate_confirmed;
			if(config.verbose>=1){
#pragma omp critical
				log("%ld\t%5ld sure closer %3d maybe closer %3d (%3d +%3d)",
						target->id,
						cand->candidates[i].mesh_wrapper->id,
						sure_closer,
						maybe_closer,
						cand->candidate_confirmed,
						cand_left);
			}
			// the rank makes sure this candidate must be a positive result
			if(maybe_closer < cand_left){
#pragma omp critical
				ctx.report_result(target->id, cand->candidates[i].mesh_wrapper->id);
				cand->candidate_confirmed++;
				//delete cand->candidates[i];
				cand->candidates.erase(cand->candidates.begin()+i);
				list_size--;
				//log("ranked %d, %d confirmed", rank, target->candidate_confirmed);
				continue;
			}

			// the rank makes sure this one should be removed as it must not be qualified
			if(sure_closer >= cand_left){
				//delete cand->candidates[i];
				cand->candidates.erase(cand->candidates.begin()+i);
				list_size--;
				continue;
			}
			i++;
		}//end for
		// the target one should be kept
	}

	// remove the queried objects and their candidates when the query result is confirmed.
	for(vector<candidate_entry *>::iterator it=ctx.candidates.begin();it!=ctx.candidates.end();){
		//cout<<(*it)->mesh_wrapper->id<<" "<<(*it)->candidate_confirmed<<" "<<ctx.knn<<endl;
		if((*it)->candidate_confirmed==config.knn){
			delete *it;
			it = ctx.candidates.erase(it);
		}else{
			it++;
		}
	}

	ctx.updatelist_time += logt("updating the candidate lists",start);
}

}
