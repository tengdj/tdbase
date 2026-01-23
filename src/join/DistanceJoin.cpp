/*
 * DistanceJoin.cpp
 *
 *  Created on: Jan 22, 2026
 *      Author: teng
 */


#include "SpatialJoin.h"

using namespace std;

namespace tdbase{


range DistanceJoin::update_voxel_pair_list(vector<voxel_pair> &voxel_pairs, double minmaxdist, bool keep_empty){

	assert(voxel_pairs.size()>0);

	if(config.verbose>=2){
		int valid_voxel = 0;
		int invalid_voxel = 0;
		for(auto &vp:voxel_pairs){
			if(vp.dist.valid()){
				valid_voxel++;
			}else{
				invalid_voxel++;
			}
		}
		log("invalid_voxles/total = %.2f\% (%d/%d)\t",invalid_voxel*1.0/(invalid_voxel+valid_voxel),invalid_voxel, (invalid_voxel+valid_voxel));
	}

	for(auto vp_iter = voxel_pairs.begin();vp_iter!=voxel_pairs.end();){
		if((vp_iter->dist.valid()&&vp_iter->dist.mindist > minmaxdist) // a closer voxel pair already exist
				||(!keep_empty&&vp_iter->has_empty_voxel())){ //remove the pairs which has an empty voxel
			// evict this unqualified voxel pairs
			voxel_pairs.erase(vp_iter);
		}else{
			vp_iter++;
		}
	}

	assert(voxel_pairs.size()>0);

	// now update the newest object-level distance from the voxel level distance
	range ret;
	ret.mindist = DBL_MAX;
	ret.maxdist = DBL_MAX;
	for(auto &vp:voxel_pairs){
		ret.mindist = min(ret.mindist, vp.dist.mindist);
		ret.maxdist = min(ret.maxdist, vp.dist.maxdist);
	}
	return ret;
}

void DistanceJoin::update_distance_ranges(query_context &ctx){

	// Called by index retrieval, need not be updated
	if(ctx.gp.results==NULL){
		return;
	}
	int index = 0;
//#pragma omp parallel for
	for(candidate_entry *c:ctx.candidates){
		HiMesh_Wrapper *wrapper1 = c->mesh_wrapper;
		for(candidate_info &ci:c->candidates){
			HiMesh_Wrapper *wrapper2 = ci.mesh_wrapper;

			if(config.use_aabb){
				// progressive query is invalid when AABB acceleration is enabled
				assert(ctx.cur_lod==config.highest_lod());
				result_container res = ctx.gp.results[index++];
				ci.distance.mindist = res.distance;
				ci.distance.maxdist = res.distance;
			}else{
				double vox_minmaxdist = DBL_MAX;
				for(voxel_pair &vp:ci.voxel_pairs){
					result_container res = ctx.gp.results[index++];
					range dist = vp.dist;

					if(vp.has_empty_voxel()){
						if(config.verbose>=1)
						{
							log("%ld(%d)\t"
								"%ld(%d):\t"
								"[%.2f, %.2f]->"
								"[%.2f, %.2f] "
								"res: [%.2f, %.2f, %.2f]",
								wrapper1->id,vp.v1->id,
								wrapper2->id,vp.v2->id,
								vp.dist.mindist, vp.dist.maxdist,
								dist.mindist, dist.maxdist,
								res.min_dist, res.distance, res.max_dist);
						}
						continue;
					}

					// update the distance
					if(ctx.cur_lod==config.highest_lod()){
						// now we have a precise distance
						dist.mindist = res.distance;
						dist.maxdist = res.distance;
					} else {
						dist.mindist = res.min_dist;
						dist.maxdist = res.max_dist;
//						dist.mindist = std::max(dist.mindist, res.min_dist);
//						dist.maxdist = std::min(dist.maxdist, res.max_dist);
					}

					if(config.verbose>=1)
					{
						log("%ld(%d)\t"
							"%ld(%d):\t"
							"[%.2f, %.2f]->"
							"[%.2f, %.2f] "
							"res: [%.2f, %.2f, %.2f]",
							wrapper1->id,vp.v1->id,
							wrapper2->id,vp.v2->id,
							vp.dist.mindist, vp.dist.maxdist,
							dist.mindist, dist.maxdist,
							res.min_dist, res.distance, res.max_dist);
					}
					vp.dist = dist;
					vox_minmaxdist = min(vox_minmaxdist, (double)dist.maxdist);
					//assert(dist.valid());
				}// end of voxel pair iteration
				// after each round, some voxels need to be evicted,
				// remove the invalid pair at the highest LOD (have an empty voxel)
				// (happens when low LOD serve as the highest LOD, some voxel will be empty)
				ci.distance = update_voxel_pair_list(ci.voxel_pairs, vox_minmaxdist,!(ctx.cur_lod==config.highest_lod()));
				assert(ci.voxel_pairs.size()>0);
				//assert(ci.distance.valid());
			}// end the if for AABB or progressive querying
		}// end evaluate the candidates for each object
	}// end evaluate all the objects' candidates
}

//utility function to calculate the distances between the facet pairs in voxel pairs in batch
void DistanceJoin::geometric_computation(query_context &ctx){
	struct timeval start = tdbase::get_cur_time();
	if(config.use_aabb){
		int index = 0;
		for(candidate_entry *c:ctx.candidates){
			for(candidate_info &info:c->candidates){
				ctx.gp.results[index++].distance = c->mesh_wrapper->get_mesh()->distance_tree(info.mesh_wrapper->get_mesh());
			}// end for distance_candiate list
		}// end for candidates
		// clear the trees for current LOD
		for(candidate_entry *c:ctx.candidates){
			for(candidate_info &info:c->candidates){
				info.mesh_wrapper->get_mesh()->clear_aabb_tree();
			}
			c->mesh_wrapper->get_mesh()->clear_aabb_tree();
		}
	}else{
		computer->get_distance(ctx.gp);
	}
	ctx.computation_time += logt("distance computation", start);
}
}
