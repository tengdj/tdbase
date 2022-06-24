/*
 * himesh_internal_index.cpp
 *
 *  Created on: Jun 24, 2022
 *      Author: teng
 */

#include "himesh.h"

namespace hispeed{

void HiMesh::get_segments(){
	segments.clear();
	for(Edge_const_iterator eit = edges_begin(); eit!=edges_end(); ++eit){
		Point p1 = eit->vertex()->point();
		Point p2 = eit->opposite()->vertex()->point();
		if(p1!=p2){
			segments.push_back(Segment(p1, p2));
		}
	}
}


void HiMesh::clear_aabb_tree(){
	if(aabb_tree){
		delete aabb_tree;
		aabb_tree = NULL;
	}
	segments.clear();
}

SegTree *HiMesh::get_aabb_tree(){

	if(aabb_tree == NULL){
		get_segments();
		aabb_tree = new SegTree(segments.begin(), segments.end());
		aabb_tree->accelerate_distance_queries();
	}
	return aabb_tree;
}

TriangleTree *HiMesh::get_aabb_tree_triangle(){
	//TriangleTree *tree = new TriangleTree(faces(*this).first, faces(*this).second, *this);
	//tree->accelerate_distance_queries();
	//return tree;
	return NULL;
}
}
