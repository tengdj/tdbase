/*
 * candidate.h
 *
 *  Created on: Sep 21, 2022
 *      Author: teng
 */

#ifndef SRC_INCLUDE_CANDIDATE_H_
#define SRC_INCLUDE_CANDIDATE_H_
#include "himesh.h"

namespace hispeed{

/*
 * one target refers to a list of candidates
 * each candidate refers to a list of candidate pairs
 * from the target and reference polyhedrons
 * */
class voxel_pair{
public:
	Voxel *v1;
	Voxel *v2;
	range dist;
	voxel_pair(Voxel *v1, Voxel *v2, range dist){
		this->v1 = v1;
		this->v2 = v2;
		this->dist = dist;
	};
	voxel_pair(Voxel *v1, Voxel *v2){
		this->v1 = v1;
		this->v2 = v2;
		this->dist = v1->distance(*v2);
	}
};

class candidate_info{
public:
	HiMesh_Wrapper *mesh_wrapper;
	range distance;
	vector<voxel_pair> voxel_pairs;
};

class candidate_entry{
public:
	candidate_entry(){}
	candidate_entry(HiMesh_Wrapper *m, vector<candidate_info> &c){
		mesh_wrapper = m;
		candidates = c;
	}
	HiMesh_Wrapper *mesh_wrapper = NULL;
	vector<candidate_info> candidates;
	int candidate_confirmed = 0;
};

size_t get_pair_num(vector<candidate_entry *> &candidates);
size_t get_candidate_num(vector<candidate_entry *> &candidates);

}


#endif /* SRC_INCLUDE_CANDIDATE_H_ */
