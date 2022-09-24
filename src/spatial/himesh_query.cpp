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
	float *tri1 = new float[9*size_of_facets()];
	fill_triangles(tri1);
	float *tri2 = new float[9*target->size_of_facets()];
	target->fill_triangles(tri2);
	float dist = TriDist_single(tri1, tri2, size_of_facets(), target->size_of_facets());
	delete []tri1;
	delete []tri2;
	return dist;
}

range HiMesh::distance_range(HiMesh *target){
	range dist;
	dist.maxdist = distance(target);
	dist.mindist = dist.maxdist - getmaximumCut() - target->getmaximumCut();
	return dist;
}

}


