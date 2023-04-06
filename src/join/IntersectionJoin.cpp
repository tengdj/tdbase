/*
 * IntersectionJoin.cpp
 *
 *  Created on: Sep 16, 2022
 *      Author: teng
 */

#include "SpatialJoin.h"

namespace hispeed{

vector<candidate_entry> SpatialJoin::mbb_intersect(Tile *tile1, Tile *tile2){
	vector<candidate_entry> candidates;
	OctreeNode *tree = tile2->build_octree(20);
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

/*
 * the main function for detecting the intersection
 * relationship among polyhedra in the tile
 *
 * */
void SpatialJoin::intersect(query_context ctx){
	struct timeval start = get_cur_time();
	struct timeval very_start = get_cur_time();

	// filtering with MBBs to get the candidate list
	vector<candidate_entry> candidates = mbb_intersect(ctx.tile1, ctx.tile2);
	ctx.index_time += hispeed::get_time_elapsed(start,false);
	logt("index retrieving", start);

	for(uint lod:ctx.lods){
		ctx.cur_lod = lod;
		struct timeval iter_start = start;
		size_t pair_num = get_pair_num(candidates);
		if(pair_num==0){
			break;
		}
		size_t candidate_num = get_candidate_num(candidates);
		log("%ld polyhedron has %d candidates %d voxel pairs %.2f voxel pairs per candidate",
				candidates.size(), candidate_num, pair_num, (1.0*pair_num)/candidates.size());

		// do the decoding, packing, and computing
		check_intersection(candidates, ctx);

		// now update the intersection status and update the all candidate list
		// report results if necessary
		int index = 0;
		start = get_cur_time();
		for(auto ce_iter=candidates.begin();ce_iter!=candidates.end();){
			HiMesh_Wrapper *wrapper1 = ce_iter->mesh_wrapper;
			//print_candidate_within(*ce_iter);
			for(auto ci_iter=ce_iter->candidates.begin();ci_iter!=ce_iter->candidates.end();){
				bool determined = false;
				HiMesh_Wrapper *wrapper2 = ci_iter->mesh_wrapper;
				for(voxel_pair &vp:ci_iter->voxel_pairs){
					determined |= ctx.results[index++].result.intersected;
				}
				//log("%d %d %d",wrapper1->id, wrapper2->id,determined);
				if(determined){
					wrapper1->report_result(wrapper2);
					ce_iter->candidates.erase(ci_iter);
				}else{
					ci_iter++;
				}
			}
			if(ce_iter->candidates.size()==0){
				candidates.erase(ce_iter);
			}else{
				ce_iter++;
			}
		}
		delete []ctx.results;
		ctx.updatelist_time += logt("update the candidate list", start);

		logt("evaluating with lod %d", iter_start, lod);
		log("");
	}
	ctx.overall_time = hispeed::get_time_elapsed(very_start, false);
//	for(auto ce_iter=candidates.begin();ce_iter!=candidates.end();ce_iter++){
//		HiMesh_Wrapper *wrapper1 = ce_iter->mesh_wrapper;
//		//print_candidate_within(*ce_iter);
//		for(auto ci_iter=ce_iter->candidates.begin();ci_iter!=ce_iter->candidates.end();ci_iter++){
//			HiMesh_Wrapper *wrapper2 = ci_iter->mesh_wrapper;
//			log("%d not intersect %d", wrapper1->id, wrapper2->id);
//		}
//	}

	for(int i=0;i<ctx.tile1->num_objects();i++){
		ctx.result_count += ctx.tile1->get_mesh_wrapper(i)->results.size();
	}
	ctx.obj_count += min(ctx.tile1->num_objects(),global_ctx.max_num_objects1);
	global_ctx.merge(ctx);
}

}

