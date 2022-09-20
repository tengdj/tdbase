/*
 * KNNJoin.cpp
 *
 *  Created on: Sep 16, 2022
 *      Author: teng
 */

#include "SpatialJoin.h"

namespace hispeed{

//#define VERBAL_3DPRO

inline float get_min_max_dist(vector<voxel_pair> &voxel_pairs){
	float minmaxdist = DBL_MAX;
	for(voxel_pair &p:voxel_pairs){
		minmaxdist = min(minmaxdist, p.dist.maxdist);
	}
	return minmaxdist;
}

inline bool update_voxel_pair_list(vector<voxel_pair> &voxel_pairs, double minmaxdist){

	int voxel_pair_size = voxel_pairs.size();
	// some voxel pair is farther than this one
	for(int j=0;j<voxel_pair_size;){
		// a closer voxel pair already exist
		if(voxel_pairs[j].dist.mindist > minmaxdist){
			// evict this dequalified voxel pairs
			voxel_pairs.erase(voxel_pairs.begin()+j);
			voxel_pair_size--;
		}else{
			j++;
		}
	}
	return true;
}

static void print_candidate(candidate_entry &cand){
#ifdef VERBAL_3DPRO
	printf("%ld (%d + %ld)\t\n", cand.first->id, cand.first->candidate_confirmed, cand.second.size());
	int i=0;
	for(candidate_info &ci:cand.second){
		printf("%d:\t%ld\t[%f,%f]\n",i++,ci.mesh_wrapper->id,ci.distance.mindist,ci.distance.maxdist);
	}
#endif
}

inline void update_candidate_list_knn(candidate_entry &cand, int knn){
	HiMesh_Wrapper *target = cand.first;
	int list_size = cand.second.size();
	for(int i=0;i<list_size && knn>target->candidate_confirmed;){
		int sure_closer = 0;
		int maybe_closer = 0;
		for(int j=0;j<cand.second.size();j++){
			if(i==j){
				continue;
			}
			// count how many candidates that are surely closer than this one
			if(cand.second[i].distance>=cand.second[j].distance) {
				sure_closer++;
			}
			// count how many candidates that are possibly closer than this one
			if(!(cand.second[i].distance<=cand.second[j].distance)) {
				maybe_closer++;
			}
		}
		int cand_left = knn-target->candidate_confirmed;
#ifdef VERBAL_3DPRO
		log("%3d %5d sure closer %3d maybe closer %3d confirmed %3d rest %3d",
				i,
				cand.second[i].mesh_wrapper->id,
				sure_closer,
				maybe_closer,
				target->candidate_confirmed,
				cand_left);
#endif
		// the rank makes sure this one is confirmed
		if(maybe_closer < cand_left){
			report_result(target->id, cand.second[i].mesh_wrapper->id);
			target->candidate_confirmed++;
			cand.second.erase(cand.second.begin()+i);
			list_size--;
			//log("ranked %d, %d confirmed", rank, target->candidate_confirmed);
			continue;
		}

		// the rank makes sure this one should be removed
		if(sure_closer >= cand_left){
			cand.second.erase(cand.second.begin()+i);
			list_size--;
			continue;
		}
		i++;
	}//end for
	// the target one should be kept
}

void evaluate_candidate_lists(vector<candidate_entry> &candidates, query_context ctx){
	for(vector<candidate_entry>::iterator it=candidates.begin();it!=candidates.end();){
		print_candidate(*it);
		update_candidate_list_knn(*it, ctx.knn);
		if(it->first->candidate_confirmed==ctx.knn){
			it = candidates.erase(it);
		}else{
			it++;
		}
	}
}

vector<candidate_entry> SpatialJoin::mbb_knn(Tile *tile1, Tile *tile2, query_context &ctx){
	vector<candidate_entry> candidates;
	vector<pair<int, range>> candidate_ids;
	OctreeNode *tree = tile2->build_octree(20);
	size_t tile1_size = min(tile1->num_objects(), ctx.max_num_objects1);
	for(int i=0;i<tile1_size;i++){
		// for each object
		//1. use the distance between the mbbs of objects as a
		//	 filter to retrieve candidate objects
		vector<candidate_info> candidate_list;
		HiMesh_Wrapper *wrapper1 = tile1->get_mesh_wrapper(i);
		float min_maxdistance = DBL_MAX;
		tree->query_knn(&(wrapper1->box), candidate_ids, min_maxdistance, ctx.knn);
		assert(!candidate_ids.size()>=ctx.knn);

		//2. we further go through the voxels in two objects to shrink
		// 	 the candidate list in a finer grain
		for(pair<int, range> &p:candidate_ids){
			HiMesh_Wrapper *wrapper2 = tile2->get_mesh_wrapper(p.first);
			candidate_info ci;
			ci.mesh_wrapper = wrapper2;
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
			update_voxel_pair_list(ci.voxel_pairs, min_maxdist);
			assert(ci.voxel_pairs.size()>0);
			// form the distance range of objects with the evaluations of voxel pairs
			ci.distance = ci.voxel_pairs[0].dist;
			for(voxel_pair &p:ci.voxel_pairs){
				ci.distance.update(p.dist);
			}
			assert(ci.distance.mindist<=ci.distance.maxdist);
			candidate_list.push_back(ci);
		}

		//log("%ld %ld", candidate_ids.size(),candidate_list.size());
		// save the candidate list
		if(candidate_list.size()>0){
			candidates.push_back(candidate_entry(wrapper1, candidate_list));
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
void SpatialJoin::nearest_neighbor(Tile *tile1, Tile *tile2, query_context ctx){
	struct timeval start = get_cur_time();
	struct timeval very_start = get_cur_time();

	// filtering with MBBs to get the candidate list
	vector<candidate_entry> candidates = mbb_knn(tile1, tile2, ctx);
	ctx.index_time += get_time_elapsed(start, false);
	logt("index retrieving", start);

	// now we start to get the distances with progressive level of details
	for(int lod:ctx.lods){
		struct timeval iter_start = get_cur_time();
		start = get_cur_time();

		// update the candidate lists with the newest distance ranges


		const int pair_num = get_pair_num(candidates);
		if(pair_num==0){
			break;
		}
		size_t candidate_num = get_candidate_num(candidates);
		log("%ld polyhedron has %d candidates %d voxel pairs %f voxel pairs per candidate",
				candidates.size(), candidate_num, pair_num, (1.0*pair_num)/candidates.size());
		// retrieve the necessary meshes
		size_t segment_pair_num = 0;

		for(candidate_entry &c:candidates){
			// the nearest neighbor is found
			assert(c.first->candidate_confirmed+c.second.size()>ctx.knn);
			print_candidate(c);
			HiMesh_Wrapper *wrapper1 = c.first;
			for(candidate_info &info:c.second){
				for(voxel_pair &vp:info.voxel_pairs){
					assert(vp.v1&&vp.v2);
					// not filled yet
					if(vp.v1->data.find(lod)==vp.v1->data.end()){
						// ensure the mesh is extracted
						tile1->decode_to(wrapper1->id, lod);
						wrapper1->fill_voxels(DT_Segment);
					}

					if(vp.v2->data.find(lod)==vp.v2->data.end()){
						tile2->decode_to(info.mesh_wrapper->id, lod);
						info.mesh_wrapper->fill_voxels(DT_Segment);
					}
					segment_pair_num += vp.v1->size[lod]*vp.v2->size[lod];
				}// end for voxel_pairs
			}// end for distance_candiate list
		}// end for candidates
		ctx.decode_time += hispeed::get_time_elapsed(start, false);
		logt("decoded with %ld segment pairs for lod %d", start, segment_pair_num, lod);
		if(segment_pair_num==0){
			log("no segments is filled in this round");
			continue;
		}
		tile1->reset_time();

		// truly conduct the geometric computations
		float *distances = calculate_distance(candidates, ctx, lod);
		ctx.computation_time += hispeed::get_time_elapsed(start, false);
		// now update the distance range with the new distances
		int index = 0;
		for(candidate_entry &c:candidates){
			for(candidate_info &ci:c.second){
				double vox_minmaxdist = DBL_MAX;
				for(voxel_pair &vp:ci.voxel_pairs){
					// update the distance
					if(vp.v1->size[lod]>0&&vp.v2->size[lod]>0){
						range dist = vp.dist;
						if(lod==ctx.highest_lod()){
							// now we have a precise distance
							dist.mindist = distances[index];
							dist.maxdist = distances[index];
						}else{
							dist.maxdist = std::min(dist.maxdist, distances[index]);
						}
						vp.dist = dist;
						if(lod==100){
							if(ci.distance.maxdist>=dist.maxdist){
								ci.distance = dist;
							}
						}else{
							ci.distance.update(dist);
						}
						vox_minmaxdist = min(vox_minmaxdist, (double)dist.maxdist);
					}
					index++;
				}
				// after each round, some voxels need to be evicted
				update_voxel_pair_list(ci.voxel_pairs, vox_minmaxdist);
			}
		}
		logt("calculating distance", start);

		// update the list after processing each LOD
		evaluate_candidate_lists(candidates, ctx);
		ctx.updatelist_time += hispeed::get_time_elapsed(start, false);
		logt("updating the candidate lists",start);

		delete distances;
		logt("evaluating with lod %d", iter_start, lod);
	}

	ctx.overall_time = hispeed::get_time_elapsed(very_start, false);
	pthread_mutex_lock(&g_lock);
	global_ctx.merge(ctx);
	pthread_mutex_unlock(&g_lock);
}

}