/*
 * himesh_query.cpp
 *
 *  Created on: Sep 24, 2022
 *      Author: teng
 */

#include "himesh.h"
#include "geometry.h"

namespace hispeed{

query_context global_ctx;

void print_triangles(float *triangle, size_t size){
	printf("OFF\n%ld %ld 0\n\n",size*3,size);
	for(size_t i=0;i<size*3;i++){
		printf("%f %f %f\n", *(triangle+3*i), *(triangle+3*i+1), *(triangle+3*i+2));
	}
	for(size_t i=0;i<size;i++){
		printf("3\t%ld %ld %ld\n",i*3,i*3+1,i*3+2);
	}
	printf("\n");
}

bool HiMesh::intersect(HiMesh *target){
	float *tri1, *tri2;
	size_t s1 = fill_triangles(tri1);
	size_t s2 = target->fill_triangles(tri2);
	result_container res = TriInt_single(tri1, tri2, s1, s2);
	if(res.result.intersected && global_ctx.verbose>=1){
		print_triangles(tri1+res.p1*9, 1);
		print_triangles(tri2+res.p2*9, 1);
	}
	delete []tri1;
	delete []tri2;
	return res.result.intersected;
}

float HiMesh::distance(HiMesh *target){
	result_container res;
	float *tri1, *tri2;
	size_t s1 = fill_triangles(tri1);
	size_t s2 = target->fill_triangles(tri2);
	res = TriDist_single(tri1, tri2, s1, s2);
	if(global_ctx.verbose>=1){
		print_triangles(tri1+res.p1*9, 1);
		print_triangles(tri2+res.p2*9, 1);
	}
	delete []tri1;
	delete []tri2;

	return res.result.distance;
}

bool HiMesh::intersect_tree(HiMesh *target){
	assert(segments.size()>0);
    for(Segment &s:segments){
    	if(target->get_aabb_tree_triangle()->do_intersect(s)){
    		return true;
    	}
    }
    return false;
}

float HiMesh::distance_tree(HiMesh *target){
	double min_dist = DBL_MAX;
	{
		list<Point> vertices = get_vertices();
		for(Point &p:vertices){
			FT sqd = target->get_aabb_tree_triangle()->squared_distance(p);
			double dist = (double)CGAL::to_double(sqd);
			if(min_dist>dist){
				min_dist = dist;
			}
		}
	}
	{
		list<Point> vertices = target->get_vertices();
		for(Point &p:vertices){
			FT sqd = get_aabb_tree_triangle()->squared_distance(p);
			double dist = (double)CGAL::to_double(sqd);
			if(min_dist>dist){
				min_dist = dist;
			}
		}
	}

	return sqrt(min_dist);
}

}


