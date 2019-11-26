/*
 * SpatialJoin.cpp
 *
 *  Created on: Nov 11, 2019
 *      Author: teng
 */


#include "SpatialJoin.h"
using namespace std;

namespace hispeed{

void SpatialJoin::formalize_computing(){

	// filtering with MBBs to get the candidate list
	vector<vector<HiMesh_Wrapper *>> candidates;
	vector<vector<range>> candidate_dist;

	struct timeval start = get_cur_time();
	int candidate_size = 0;
	for(int i=0;i<tile1->num_objects();i++){
		vector<HiMesh_Wrapper *> geom_list;
		vector<range> range_list;

		aab b1 = tile1->get_mbb(i);
		for(int j=0;j<tile2->num_objects();j++){
			if(tile1==tile2&&i==j){
				continue;
			}
			HiMesh_Wrapper *wrapper = tile2->get_mesh_wrapper(j);
			range d = b1.distance(wrapper->box);
			// now update the list
			bool keep = true;
			int list_size = range_list.size();
			for(int i=0;i<list_size;){
				// should not keep since there is a closer one
				if(d>range_list[i]){
					keep = false;
					break;
				// one in the list cannot be the closest because of this one
				}else if(d<range_list[i]){
					range_list.erase(range_list.begin()+i);
					geom_list.erase(geom_list.begin()+i);
					list_size--;
				}else{
					i++;
				}
			}
			if(keep){
				range_list.push_back(d);
				geom_list.push_back(wrapper);
			}
		}
		candidate_size += geom_list.size();
		// save the candidate list
		candidates.push_back(geom_list);
		candidate_dist.push_back(range_list);
	}
	fprintf(stderr, "%d polyhedron has %d candidates takes %f ms\n",
			tile1->num_objects(), candidate_size, get_time_elapsed(start, true));

	int lod = 20;
	// retrieve the meshes first
	long total_segments = 0;
	for(int i=0;i<tile1->num_objects();i++){
		HiMesh *mesh1 = tile1->get_mesh(i, lod);
		for(int j=0;j<candidates[i].size();j++){
			HiMesh *mesh2 = tile2->get_mesh(candidates[i][j]->id, lod);
			total_segments += mesh1->get_segment_num();
			total_segments += mesh2->get_segment_num();
		}
	}
	report_time("getting meshes", start);

	// organize the data for computing
	float *data = new float[6*total_segments];
	long *offset_size = new long[3*candidate_size];
	float *distances = new float[candidate_size];
	memset(distances, 0, candidate_size*sizeof(float));
	long offset = 0;
	int index = 0;
	for(int i=0;i<tile1->num_objects();i++){
		HiMesh *mesh1 = tile1->get_mesh(i, lod);
		for(int j=0;j<candidates[i].size();j++){
			HiMesh *mesh2 = tile2->get_mesh(candidates[i][j]->id, lod);
			offset_size[3*index] = offset;
			mesh1->get_segments(&data[offset]);
			offset_size[3*index+1] = mesh1->get_segment_num();
			offset += offset_size[3*index+1];
			mesh2->get_segments(&data[offset]);
			offset_size[3*index+2] = mesh2->get_segment_num();
			offset += offset_size[3*index+2];
			index++;
		}
	}
	assert(index==candidate_size);
	report_time("organizing data", start);
	hispeed::SegDist_batch_gpu(data, offset_size, distances, candidate_size);
	report_time("get distance with GPU", start);
	index = 0;
	for(int i=0;i<tile1->num_objects();i++){
		for(int j=0;j<candidates[i].size();j++){
			cout<<candidate_dist[i][j]<<" "<<distances[index++]<<endl;
		}
		break;
	}
	delete data;
	delete offset_size;
	delete distances;
}


}


