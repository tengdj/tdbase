/*
 * test.cpp
 *
 *  Created on: Oct 24, 2019
 *      Author: teng
 */

#include <stdio.h>
#include <iostream>
#include <type_traits>
#include <math.h>
#include "util/util.h"
using namespace std;
using namespace hispeed;

inline void normalize_vec(float vec[3]){
	float v = 0;
	for(int i=0;i<3;i++){
		v += vec[i]*vec[i];
	}
	v = sqrt(v);
	for(int i=0;i<3;i++){
		vec[i] = vec[i]/v;
	}
}

inline void cross_product(float vec1[3], float vec2[3], float result[3]){
	result[0] = vec1[1]*vec2[2]-vec1[2]*vec2[1];
	result[1] = vec1[2]*vec2[0]-vec1[0]*vec2[2];
	result[2] = vec1[0]*vec2[1]-vec1[1]*vec2[0];
}

int main(int argc, char **argv){
	float vec1[3],vec2[3],result[3];
	vec1[0] = -1.0;
	vec1[1] = 0;
	vec1[2] = 1.0;
	vec2[0] = 1.0;
	vec2[1] = 0.0;
	vec2[2] = 1.0;
	normalize_vec(vec1);
	normalize_vec(vec2);

	cross_product(vec1,vec2,result);
	cout<<result[0]<<" "<<result[1]<<" "<<result[2]<<endl;
	for(int i=0;i<3;i++){
		vec1[i] = -vec1[i];
	}
	cross_product(vec1,vec2,result);
	cout<<result[0]<<" "<<result[1]<<" "<<result[2]<<endl;
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
