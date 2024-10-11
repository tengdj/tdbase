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

}
