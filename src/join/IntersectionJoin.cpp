/*
 * IntersectionJoin.cpp
 *
 *  Created on: Sep 16, 2022
 *      Author: teng
 */

#include "SpatialJoin.h"

namespace tdbase{

void IntersectJoin::index_retrieval(Tile *tile1, Tile *tile2, query_context &ctx){
	struct timeval start = get_cur_time();

	OctreeNode *tree = tile2->get_octree();
#pragma omp parallel for
	for(int i=0;i<tile1->num_objects();i++){
		vector<int> candidate_ids;
		HiMesh_Wrapper *wrapper1 = tile1->get_mesh_wrapper(i);

		tree->query_intersect(&(wrapper1->box), candidate_ids);
		// no candidates
		if(candidate_ids.empty()){
			continue;
		}
		candidate_entry *ce = new candidate_entry(wrapper1);

		std::sort(candidate_ids.begin(), candidate_ids.end());
		int former = -1;
		for(int tile2_id:candidate_ids){
			if(tile2_id==former){
				// duplicate
				continue;
			}
			HiMesh_Wrapper *wrapper2 = tile2->get_mesh_wrapper(tile2_id);
			candidate_info ci(wrapper2);
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
				ce->add_candidate(ci);
			}else{
				//delete ci;
			}
			former = tile2_id;
		}
		candidate_ids.clear();
		// save the candidate list if needed
		if(ce->candidates.size()>0){
#pragma omp critical
			ctx.candidates.push_back(ce);
		}else{
			delete ce;
		}
		candidate_ids.clear();
	}
	ctx.index_time += logt("index retrieving", start);
}

void IntersectJoin::geometric_computation(query_context &ctx){
	struct timeval start = tdbase::get_cur_time();
	if(config.use_aabb){
		int index = 0;
		for(candidate_entry *c:ctx.candidates){
			c->mesh_wrapper->get_mesh()->get_segments();
			for(candidate_info &info:c->candidates){
				assert(info.voxel_pairs.size()==1);
				ctx.gp.results[index++].intersected = c->mesh_wrapper->get_mesh()->intersect_tree(info.mesh_wrapper->get_mesh());
			}// end for candidate list
		}// end for candidates
		// clear the trees for current LOD
		for(candidate_entry *c:ctx.candidates){
			for(candidate_info &info:c->candidates){
				info.mesh_wrapper->get_mesh()->clear_aabb_tree();
			}
		}
	}else{
		computer->get_intersect(ctx.gp);
	}
	ctx.computation_time += logt("computation for checking intersection", start);
}

void IntersectJoin::evaluate_candidate_lists(query_context &ctx){
	struct timeval start = get_cur_time();

	int index = 0;
	// update the candidates with the calculated intersection info
	for(auto ce_iter=ctx.candidates.begin();ce_iter!=ctx.candidates.end();){
		HiMesh_Wrapper *wrapper1 = (*ce_iter)->mesh_wrapper;
		//print_candidate_within(*ce_iter);
		for(auto ci_iter=(*ce_iter)->candidates.begin();ci_iter!=(*ce_iter)->candidates.end();){
			HiMesh_Wrapper *wrapper2 = (ci_iter)->mesh_wrapper;
			bool determined = false;
			int cand_count = 0;
			for(voxel_pair &vp:(ci_iter)->voxel_pairs){
				determined |= ctx.gp.results[index].intersected;
				// the minimum possible distance already been computed
				cand_count += (ctx.gp.results[index].min_dist>0);
				index++;
			}

			if(determined){
				// must intersect
				ctx.report_result(wrapper1->id, wrapper2->id);

				//delete *ci_iter;
				(*ce_iter)->candidates.erase(ci_iter);
				// all voxel pairs must not intersect
			}else if(cand_count == (ci_iter)->voxel_pairs.size()){
				// must not intersect
				(*ce_iter)->candidates.erase(ci_iter);
			}else{
				ci_iter++;
			}
		}
		// relationship with all candidates are determined
		if((*ce_iter)->candidates.size()==0){
			delete (*ce_iter);
			ctx.candidates.erase(ce_iter);
		}else{
			ce_iter++;
		}
	}
	ctx.updatelist_time += logt("updating the candidate lists",start);

}

/////*
//// * the main function for detecting the intersection
//// * relationship among polyhedra in the tile
//// *
//// * */
//void SpatialJoin::intersect(Tile *tile1, Tile *tile2){
//	struct timeval start = get_cur_time();
//	struct timeval very_start = get_cur_time();
//
//	query_context ctx;
//
//	// filtering with MBBs to get the candidate list
//	vector<candidate_entry *> candidates = mbb_intersect(tile1, tile2);
//	ctx.index_time += tdbase::get_time_elapsed(start,false);
//	logt("index retrieving", start);
//
//	for(uint32_t lod:config.lods){
//		ctx.cur_lod = lod;
//		struct timeval iter_start = start;
//		size_t pair_num = get_pair_num(candidates);
//		if(pair_num==0){
//			break;
//		}
//		size_t candidate_num = get_candidate_num(candidates);
//		log("%ld polyhedron has %d candidates %d voxel pairs %.2f voxel pairs per candidate",
//				candidates.size(), candidate_num, pair_num, (1.0*pair_num)/candidates.size());
//
//		// do the decoding, packing, and computing
//		check_intersection(candidates, ctx);
//
//		// now update the intersection status and update the all candidate list
//		// report results if necessary
//		start = get_cur_time();
//
//		int index = 0;
//		// update the candidates with the calculated intersection info
//		for(auto ce_iter=candidates.begin();ce_iter!=candidates.end();){
//			HiMesh_Wrapper *wrapper1 = (*ce_iter)->mesh_wrapper;
//			//print_candidate_within(*ce_iter);
//			for(auto ci_iter=(*ce_iter)->candidates.begin();ci_iter!=(*ce_iter)->candidates.end();){
//				bool determined = false;
//				HiMesh_Wrapper *wrapper2 = (ci_iter)->mesh_wrapper;
//				int cand_count = 0;
//				for(voxel_pair &vp:(ci_iter)->voxel_pairs){
//					determined |= ctx.gp.results[index].intersected;
//					// the minimum possible distance already been computed
//					cand_count += (ctx.gp.results[index].min_dist>0);
//					index++;
//				}
//
//				//log("%d %d %d",wrapper1->id, wrapper2->id,determined);
//				if(determined){
//					// must intersect
//					ctx.report_result(wrapper1->id, wrapper2->id);
//
//					//delete *ci_iter;
//					(*ce_iter)->candidates.erase(ci_iter);
//					// all voxel pairs must not intersect
//				}else if(cand_count == (ci_iter)->voxel_pairs.size()){
//					// must not intersect
//					(*ce_iter)->candidates.erase(ci_iter);
//				}else{
//					ci_iter++;
//				}
//			}
//			if((*ce_iter)->candidates.size()==0){
//				delete (*ce_iter);
//				candidates.erase(ce_iter);
//			}else{
//				ce_iter++;
//			}
//		}
//		delete []ctx.gp.results;
//		ctx.updatelist_time += logt("update the candidate list", start);
//
//		logt("evaluating with lod %d", iter_start, lod);
//		log("");
//	}
//	ctx.overall_time = tdbase::get_time_elapsed(very_start, false);
//	ctx.obj_count = tile1->num_objects();
//	if(config.print_result){
//		ctx.print_result();
//	}
//	ctx.report();
//}

}

