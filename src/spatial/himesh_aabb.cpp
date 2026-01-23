/*
 * himesh_internal_index.cpp
 *
 *  Created on: Jun 24, 2022
 *      Author: teng
 */

#include "himesh.h"

namespace tdbase{

void HiMesh::clear_aabb_tree(){
	if(triangle_tree){
		delete triangle_tree;
		triangle_tree = NULL;
	}
	aabb_triangles.clear();
}

TriangleTree *HiMesh::get_aabb_tree_triangle(){
	if(triangle_tree == NULL){
		updateAABB();
	}
	return triangle_tree;
}

void HiMesh::updateAABB(){
	if(triangle_tree){
		aabb_triangles.clear();
		delete triangle_tree;
		triangle_tree = NULL;
	}
	aabb_triangles = get_triangles();
	triangle_tree = new TriangleTree(aabb_triangles.begin(), aabb_triangles.end());
	triangle_tree->build();
	triangle_tree->accelerate_distance_queries();
}

}
