/*
 * KNNJoin.cpp
 *
 *  Created on: Sep 16, 2022
 *      Author: teng
 */

#include "SpatialJoin.h"

namespace tdbase{

inline float get_min_max_dist(vector<voxel_pair> &voxel_pairs){
	float minmaxdist = DBL_MAX;
	for(voxel_pair &p:voxel_pairs){
		minmaxdist = min(minmaxdist, p.dist.maxdist);
	}
	return minmaxdist;
}

inline range update_voxel_pair_list(vector<voxel_pair> &voxel_pairs, double minmaxdist){

	range ret;
	ret.mindist = DBL_MAX;
	ret.maxdist = minmaxdist;
	// some voxel pair is farther than this one
	for(auto vp_iter = voxel_pairs.begin();vp_iter!=voxel_pairs.end();){
		// a closer voxel pair already exist
		if(vp_iter->dist.mindist > minmaxdist){
			// evict this unqualified voxel pairs
			voxel_pairs.erase(vp_iter);
		}else{
			ret.mindist = min(ret.mindist, vp_iter->dist.mindist);
			vp_iter++;
		}
	}
	return ret;
}

void print_candidate(candidate_entry *cand){
	if(global_ctx.verbose>=1){
		log("%ld (%d + %ld)", cand->mesh_wrapper->id, cand->candidate_confirmed, cand->candidates.size());
		int i=0;
		for(candidate_info &ci:cand->candidates){
			log("%d\t%ld:\t[%f,%f]",i++,ci.mesh_wrapper->id,ci.distance.mindist,ci.distance.maxdist);
		}
	}
}

inline void update_candidate_list_knn(candidate_entry *cand, query_context &ctx){
	HiMesh_Wrapper *target = cand->mesh_wrapper;
	int list_size = cand->candidates.size();

	for(int i=0;i<list_size && ctx.knn>cand->candidate_confirmed;){
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
		int cand_left = ctx.knn-cand->candidate_confirmed;
		if(global_ctx.verbose>=1){
			log("%ld\t%5ld sure closer %3d maybe closer %3d (%3d +%3d)",
					cand->mesh_wrapper->id,
					cand->candidates[i].mesh_wrapper->id,
					sure_closer,
					maybe_closer,
					cand->candidate_confirmed,
					cand_left);
		}
		// the rank makes sure this one is confirmed
		if(maybe_closer < cand_left){
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

bool result_sort(pair<int, int> a, pair<int, int> b){
	if(a.first<b.first){
		return true;
	}else if(a.first>b.first){
		return false;
	}else{
		return a.second<=b.second;
	}
}

void evaluate_candidate_lists(vector<candidate_entry *> &candidates, query_context &ctx){

//#pragma omp parallel for
	for (int i = 0; i < candidates.size();i++) {
		update_candidate_list_knn(candidates[i], ctx);
	}

	// query result confirmed.
	for(vector<candidate_entry *>::iterator it=candidates.begin();it!=candidates.end();){
		//cout<<(*it)->mesh_wrapper->id<<" "<<(*it)->candidate_confirmed<<" "<<ctx.knn<<endl;
		if((*it)->candidate_confirmed==ctx.knn){
			delete *it;
			it = candidates.erase(it);
		}else{
			it++;
		}
	}
}

vector<candidate_entry *> SpatialJoin::mbb_knn(Tile *tile1, Tile *tile2, query_context &ctx){
	vector<candidate_entry *> candidates;
	OctreeNode *tree = tile2->get_octree();
	size_t tile1_size = min(tile1->num_objects(), ctx.max_num_objects1);

//#pragma omp parallel for
	for(int i=0;i<tile1_size;i++){
		vector<pair<int, range>> candidate_ids;
		// for each object
		//1. use the distance between the mbbs of objects as a
		//	 filter to retrieve candidate objects
		HiMesh_Wrapper *wrapper1 = tile1->get_mesh_wrapper(i);
		float min_maxdistance = DBL_MAX;
		tree->query_knn(&(wrapper1->box), candidate_ids, min_maxdistance, ctx.knn);
		assert(candidate_ids.size()>=ctx.knn);

		// result determined with only evaluating the MBBs
		if(candidate_ids.size() == ctx.knn){
			for(pair<int, range> &p:candidate_ids){
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
			ci.distance = update_voxel_pair_list(ci.voxel_pairs, min_maxdist);
			assert(ci.voxel_pairs.size()>0);
			assert(ci.distance.valid());
			ce->add_candidate(ci);
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

	// the candidates list need be evaluated after checking with the mbb
	// some queries might be answered with only querying the index
	evaluate_candidate_lists(candidates, ctx);

	return candidates;
}

/*
 * the main function for getting the nearest neighbor
 *
 * */

void SpatialJoin::nearest_neighbor(query_context ctx){
	struct timeval start = get_cur_time();
	struct timeval very_start = get_cur_time();

	// filtering with MBBs to get the candidate list
	vector<candidate_entry *> candidates = mbb_knn(ctx.tile1, ctx.tile2, ctx);
	ctx.index_time += logt("index retrieving", start);

	// now we start to get the distances with progressive level of details
	for(uint32_t lod:ctx.lods){

		ctx.cur_lod = lod;
		struct timeval iter_start = get_cur_time();
		start = get_cur_time();

		const int pair_num = get_pair_num(candidates);
		if(pair_num==0){
			break;
		}
		size_t candidate_num = get_candidate_num(candidates);
		log("%ld polyhedron has %d candidates %d voxel pairs %.2f voxel pairs per candidate",
				candidates.size(), candidate_num, pair_num, (1.0*pair_num)/candidates.size());

		// truly conduct the geometric computations
		calculate_distance(candidates, ctx);

		// now update the distance range with the new distances
		int index = 0;
		start = get_cur_time();
//#pragma omp parallel for
		for(candidate_entry *c:candidates){
			HiMesh_Wrapper *wrapper1 = c->mesh_wrapper;
			for(candidate_info &ci:c->candidates){
				HiMesh_Wrapper *wrapper2 = ci.mesh_wrapper;

				if(ctx.use_aabb){
					// progressive query is invalid when AABB acceleration is enabled
					assert(lod==ctx.highest_lod());
					result_container res = ctx.tmp_results[index++];
					ci.distance.mindist = res.distance;
					ci.distance.maxdist = res.distance;
				}else{
					double vox_minmaxdist = DBL_MAX;
					for(voxel_pair &vp:ci.voxel_pairs){
						result_container res = ctx.tmp_results[index++];
						// update the distance
						if(!vp.has_empty_voxel()){
							range dist = vp.dist;
							if(lod==ctx.highest_lod()){
//								if(!dist.valid()){
//									printf("%f\t%f\t%f\n",dist.mindist,dist.maxdist,res.distance);
//								}
								// now we have a precise distance
								dist.mindist = res.distance;
								dist.maxdist = res.distance;
							} else {
								dist.mindist = std::max(dist.mindist, res.min_dist);
								dist.maxdist = std::min(dist.maxdist, res.max_dist);
							}
							//dist.maxdist = std::max(dist.maxdist, dist.mindist);

							if(global_ctx.verbose>=1&&!dist.valid())
							{
								pid_t tid = gettid();
								log("%ld\t"
										"%ld(%d)\t"
										"%ld(%d):\t"
										"[%.2f, %.2f]->"
										"[%.2f, %.2f] "
										"res: [%.2f, %.2f, %.2f]",
										tid,
										wrapper1->id,res.p1,
										wrapper2->id,res.p2,
										vp.dist.mindist, vp.dist.maxdist,
										dist.mindist, dist.maxdist,
										res.min_dist, res.distance, res.max_dist);
							}
							vp.dist = dist;
							vox_minmaxdist = min(vox_minmaxdist, (double)dist.maxdist);
							//assert(dist.valid());
						}
					}// end of voxel pair iteration
					// after each round, some voxels need to be evicted,
					// remove the invalid pair at the highest LOD (have an empty voxel)
					// (happens when low LOD serve as the highest LOD, some voxel will be empty)
					ci.distance = update_voxel_pair_list(ci.voxel_pairs, vox_minmaxdist,lod!=ctx.highest_lod());
					assert(ci.voxel_pairs.size()>0);
					assert(ci.distance.mindist<=ci.distance.maxdist);
				}// end the if for AABB or progressive querying
			}// end evaluate the candidates for each object
		}// end evaluate all the objects' candidates

		// update the list after processing each LOD
		evaluate_candidate_lists(candidates, ctx);
		delete []ctx.tmp_results;
		ctx.updatelist_time += logt("updating the candidate lists",start);

		logt("evaluating with lod %d. %ld candidates left", iter_start, lod,candidates.size());
		log("");
		if(lod==ctx.highest_lod()){
			assert(candidates.size()==0);
		}
	}

	ctx.overall_time = tdbase::get_time_elapsed(very_start, false);
	ctx.obj_count += min(ctx.tile1->num_objects(),global_ctx.max_num_objects1);
	global_ctx.merge(ctx);
}

}
