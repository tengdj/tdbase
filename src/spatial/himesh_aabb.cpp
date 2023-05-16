/*
 * himesh_internal_index.cpp
 *
 *  Created on: Jun 24, 2022
 *      Author: teng
 */

#include "himesh.h"

namespace hispeed{

void HiMesh::clear_aabb_tree(){
	if(segment_tree){
		delete segment_tree;
		segment_tree = NULL;
	}
	segments.clear();
	if(triangle_tree){
		delete triangle_tree;
		triangle_tree = NULL;
	}
	triangles.clear();
}

SegTree *HiMesh::get_aabb_tree_segment(){
	if(segment_tree == NULL){
		get_segments();
		segment_tree = new SegTree(segments.begin(), segments.end());
		segment_tree->build();
		segment_tree->accelerate_distance_queries();
	}
	return segment_tree;
}

TriangleTree *HiMesh::get_aabb_tree_triangle(){
	if(triangle_tree == NULL){
		//struct timeval start = get_cur_time();
		get_triangles();
		triangle_tree = new TriangleTree(triangles.begin(), triangles.end());
		triangle_tree->build();
		triangle_tree->accelerate_distance_queries();
		//logt("building triangle tree", start);
	}
	return triangle_tree;
}

//Tree *HiMesh::get_aabb_tree(){
//    Tree tree(faces(*this).first, faces(*this).second, *this);
//    // query point
//    Point query(0.0, 0.0, 3.0);
//    // computes squared distance from query
//    FT sqd = tree.squared_distance(query);
//    std::cout << "squared distance: " << sqd << std::endl;
//    // computes closest point
//    Point closest = tree.closest_point(query);
//    std::cout << "closest point: " << closest << std::endl;
//    // computes closest point and primitive id
//    Point_and_primitive_id pp = tree.closest_point_and_primitive(query);
//    Point closest_point = pp.first;
//    Polyhedron::Face_handle f = pp.second; // closest primitive id
//
//    std::cout << "closest point: " << closest_point << std::endl;
//    std::cout << "closest triangle: ( "
//              << f->halfedge()->vertex()->point() << " , "
//              << f->halfedge()->next()->vertex()->point() << " , "
//              << f->halfedge()->next()->next()->vertex()->point()
//              << " )" << std::endl;
//
//
//    return NULL;
//}
}
