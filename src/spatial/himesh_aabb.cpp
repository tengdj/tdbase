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


void HiMesh::get_triangles(){
	triangles.clear();
	size_t size = size_of_facets();
	for ( Facet_const_iterator f = facets_begin(); f != facets_end(); ++f){
		Point p1 = f->halfedge()->vertex()->point();
		Point p2 = f->halfedge()->next()->vertex()->point();
		Point p3 = f->halfedge()->next()->next()->vertex()->point();
		triangles.push_back(Triangle(p1, p2, p3));
	}
}

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
		get_triangles();
		triangle_tree = new TriangleTree(triangles.begin(), triangles.end());
		triangle_tree->build();
		triangle_tree->accelerate_distance_queries();
	}
	return triangle_tree;
}


}
