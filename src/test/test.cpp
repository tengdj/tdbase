/*
 * test.cpp
 *
 *  Created on: Oct 24, 2019
 *      Author: teng
 */

#include <stdio.h>
#include <iostream>
#include <type_traits>
#include "util/util.h"
using namespace std;
using namespace hispeed;
int main(int argc, char **argv){
	struct timeval t = get_cur_time();
	int i = 10;
	logt("i=%d",t,i);
}



//if(false){
//	HiMesh *mesh = tile2->get_mesh(0,100);
//	HiMesh *mesh2 = tile2->get_mesh(1,100);
//	logt("get meshs", start);
//	SegTree *tree = mesh->get_aabb_tree();
//	logt("get aabb tree", start);
//	std::vector<Point> points;
//	mesh2->get_vertices(points);
//	float mindist = DBL_MAX;
//	for(Point p:points){
//		float dist = (float)CGAL::to_double(tree->squared_distance(p));
//		if(dist<mindist){
//			mindist = dist;
//		}
//	}
//	cout<<sqrt(mindist)<<endl;
//	logt("getting min", start);
//}
//
//if(false){
//	TriangleTree *tree1 = get_aabb_tree(hispeed::read_polyhedron());
//	TriangleTree *tree2 = get_aabb_tree(hispeed::read_polyhedron());
//	logt("get aabb tree", start);
//	for(int i=0;i<tile1->num_objects();i++){
//		float mindist1 = DBL_MAX;
//		float mindist2 = DBL_MAX;
//		HiMesh *mesh = tile1->get_mesh(i,100);
//		std::vector<Point> points;
//		mesh->get_vertices(points);
//		for(Point p:points){
//			float dist = (float)CGAL::to_double(tree1->squared_distance(p));
//			if(dist<mindist1){
//				mindist1 = dist;
//			}
//			dist = (float)CGAL::to_double(tree2->squared_distance(p));
//			if(dist<mindist2){
//				mindist2 = dist;
//			}
//		}
//		points.clear();
//	}
//	logt("getting min", start);
//}
