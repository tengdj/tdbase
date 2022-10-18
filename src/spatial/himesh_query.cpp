/*
 * himesh_query.cpp
 *
 *  Created on: Sep 24, 2022
 *      Author: teng
 */

#include "himesh.h"
#include "geometry.h"

namespace hispeed{

bool HiMesh::intersect(HiMesh *target){
	float *tri1 = new float[9*size_of_facets()];
	fill_triangles(tri1);
	float *tri2 = new float[9*target->size_of_facets()];
	target->fill_triangles(tri2);
	bool inter = TriInt_single(tri1, tri2, size_of_facets(), target->size_of_facets());
	delete []tri1;
	delete []tri2;
	return inter;
}

float HiMesh::distance(HiMesh *target){
	float dist;
	if(global_ctx.etype == DT_Triangle){
		float *tri1 = new float[9*size_of_facets()];
		assert(fill_triangles(tri1) == size_of_facets());
		float *tri2 = new float[9*target->size_of_facets()];
		assert(target->fill_triangles(tri2) == target->size_of_facets());
		dist = TriDist_single(tri1, tri2, size_of_facets(), target->size_of_facets());

//		for(size_t i=0;i<size_of_facets();i++){
//			printf("%f %f %f\n", *(tri1+3*i), *(tri1+3*i+1), *(tri1+3*i+2));
//		}
//		printf("\n");
//		for(size_t i=0;i<target->size_of_facets();i++){
//			printf("%f %f %f\n", *(tri2+3*i), *(tri2+3*i+1), *(tri2+3*i+2));
//		}
		delete []tri1;
		delete []tri2;
	}else{
		float *seg1 = new float[6*size_of_edges()];
		assert(fill_segments(seg1) == size_of_edges());
		float *seg2 = new float[6*target->size_of_edges()];
		assert(target->fill_segments(seg2) == target->size_of_edges());
		dist = SegDist_single(seg1, seg2, size_of_edges(), target->size_of_edges());
		delete []seg1;
		delete []seg2;
	}
	return dist;
}

float HiMesh::distance_tree(HiMesh *target){
	double min_dist = DBL_MAX;
	{
		vector<Point> vertices;
		get_vertices(vertices);
		for(Point &p:vertices){
			FT sqd;
			if(global_ctx.etype == DT_Segment){
				sqd = target->get_aabb_tree_segment()->squared_distance(p);
			}else{
				sqd = target->get_aabb_tree_triangle()->squared_distance(p);
			}
			double dist = (double)CGAL::to_double(sqd);
			if(min_dist>dist){
				min_dist = dist;
			}
		}
		vertices.clear();
	}
	{
		vector<Point> vertices;
		target->get_vertices(vertices);
		for(Point &p:vertices){
			FT sqd;
			if(global_ctx.etype == DT_Segment){
				sqd = get_aabb_tree_segment()->squared_distance(p);
			}else{
				sqd = get_aabb_tree_triangle()->squared_distance(p);
			}
			double dist = (double)CGAL::to_double(sqd);
			if(min_dist>dist){
				min_dist = dist;
			}
		}
		vertices.clear();
	}

	return sqrt(min_dist);
}

range HiMesh::distance_range(HiMesh *target){
	range dist;
	dist.maxdist = distance(target);
	dist.mindist = dist.maxdist - getmaximumCut() - target->getmaximumCut();
	return dist;
}

}


