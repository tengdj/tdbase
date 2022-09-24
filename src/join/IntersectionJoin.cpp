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
void SpatialJoin::intersect(Tile *tile1, Tile *tile2, query_context ctx){
	struct timeval start = get_cur_time();
	struct timeval very_start = get_cur_time();

	// filtering with MBBs to get the candidate list
	vector<candidate_entry> candidates = mbb_intersect(tile1, tile2);
	ctx.index_time += hispeed::get_time_elapsed(start,false);
	logt("comparing mbbs", start);

	// now we start to ensure the intersection with progressive level of details
	size_t triangle_pair_num = 0;

	for(int lod:ctx.lods){
		struct timeval iter_start = start;
		size_t pair_num = get_pair_num(candidates);
		if(pair_num==0){
			break;
		}
		log("%ld polyhedron has %ld candidates", candidates.size(), pair_num);
		// retrieve the necessary meshes
		map<Voxel *, std::pair<uint, uint>> voxel_map;
		size_t triangle_num = 0;

		for(candidate_entry &c:candidates){
			HiMesh_Wrapper *wrapper1 = c.mesh_wrapper;
			for(candidate_info &info:c.candidates){
				HiMesh_Wrapper *wrapper2 = info.mesh_wrapper;
				for(voxel_pair vp:info.voxel_pairs){
					// not filled yet
					if(vp.v1->data.find(lod)==vp.v1->data.end())
					{
						tile1->decode_to(wrapper1->id, lod);
						wrapper1->fill_voxels(DT_Triangle);
					}
					if(vp.v2->data.find(lod)==vp.v2->data.end())
					{
						tile2->decode_to(wrapper2->id, lod);
						wrapper2->fill_voxels(DT_Triangle);
					}

					// update the voxel map
					for(int i=0;i<2;i++){
						Voxel *tv = i==0?vp.v1:vp.v2;
						if(voxel_map.find(tv)==voxel_map.end()){
							std::pair<uint, uint> p;
							p.first = triangle_num;
							p.second = tv->size[lod];
							triangle_num += tv->size[lod];
							voxel_map[tv] = p;
						}
					}

					triangle_pair_num += vp.v1->size[lod]*vp.v2->size[lod];
				}// end for voxel_pairs
			}// end for distance_candiate list
		}// end for candidates
		ctx.decode_time += hispeed::get_time_elapsed(start, false);
		logt("decoded %ld voxels with %ld triangles %ld pairs for lod %d",
				start, voxel_map.size(), triangle_num, triangle_pair_num, lod);

		tile1->reset_time();
		tile2->reset_time();
		// now we allocate the space and store the data in a buffer
		float *data = new float[9*triangle_num];
		for (map<Voxel *, std::pair<uint, uint>>::iterator it=voxel_map.begin();
				it!=voxel_map.end(); ++it){
			if(it->first->size[lod]>0){
				memcpy(data+it->second.first*9, it->first->data[lod], it->first->size[lod]*9*sizeof(float));
			}
		}
		// organize the data for computing
		uint *offset_size = new uint[4*pair_num];
		uint *intersect_status = new uint[pair_num];
		for(int i=0;i<pair_num;i++){
			intersect_status[i] = 0;
		}
		int index = 0;
		for(candidate_entry c:candidates){
			for(candidate_info &info:c.candidates){
				for(voxel_pair &vp:info.voxel_pairs){
					offset_size[4*index] = voxel_map[vp.v1].first;
					offset_size[4*index+1] = voxel_map[vp.v1].second;
					offset_size[4*index+2] = voxel_map[vp.v2].first;
					offset_size[4*index+3] = voxel_map[vp.v2].second;
					index++;
				}
			}
		}
		assert(index==pair_num);
		ctx.packing_time += hispeed::get_time_elapsed(start, false);
		logt("organizing data", start);
		geometry_param gp;
		gp.data = data;
		gp.pair_num = pair_num;
		gp.offset_size = offset_size;
		gp.intersect = intersect_status;
		gp.data_size = triangle_num;
		computer->get_intersect(gp);
		ctx.computation_time += hispeed::get_time_elapsed(start, false);
		logt("checking intersection", start);

		// now update the intersection status and update the all candidate list
		// report results if necessary
		index = 0;
		for(auto ce_iter=candidates.begin();ce_iter!=candidates.end();){
			HiMesh_Wrapper *wrapper1 = ce_iter->mesh_wrapper;
			//print_candidate_within(*ce_iter);
			for(auto ci_iter=ce_iter->candidates.begin();ci_iter!=ce_iter->candidates.end();){
				bool determined = false;
				HiMesh_Wrapper *wrapper2 = ci_iter->mesh_wrapper;
				for(voxel_pair &vp:ci_iter->voxel_pairs){
					determined |= intersect_status[index++];
				}
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

		delete []data;
		delete []offset_size;
		delete []intersect_status;
		voxel_map.clear();
		logt("current iteration", iter_start);
	}
	ctx.overall_time = hispeed::get_time_elapsed(very_start, false);

	for(int i=0;i<tile1->num_objects();i++){
		ctx.result_count += tile1->get_mesh_wrapper(i)->results.size();
	}
	ctx.obj_count += min(tile1->num_objects(),global_ctx.max_num_objects1);
	global_ctx.merge(ctx);
}

}

