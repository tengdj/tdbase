/*
 * getoff.cpp
 *
 *  Created on: Nov 20, 2019
 *      Author: teng
 */

#include "../spatial/spatial.h"
#include "../spatial/himesh.h"

using namespace hispeed;

Polyhedron adjust_polyhedron(float shift[3], float shrink, Polyhedron &poly_o){

	Polyhedron poly;
	stringstream ss;
	ss << poly_o;
	ss >> poly;

	double min_co[3] = {DBL_MAX,DBL_MAX,DBL_MAX};
	double max_co[3] = {0,0,0};

	for(Polyhedron::Vertex_iterator vi=poly.vertices_begin();vi!=poly.vertices_end();vi++){
		Point p = vi->point();
		for(int i=0;i<3;i++){
			if(min_co[i]>p[i]){
				min_co[i] = p[i];
			}
			if(max_co[i]<p[i]){
				max_co[i] = p[i];
			}
		}
	}

	for(Polyhedron::Vertex_iterator vi=poly.vertices_begin();vi!=poly.vertices_end();vi++){
		Point p = vi->point();
		if(shrink>1){
			vi->point() = Point((p[0]+shift[0]+min_co[0])/shrink, (p[1]+shift[1]+min_co[1])/shrink, (p[2]+shift[2]+min_co[2])/shrink);
		}else{
			vi->point() = Point((p[0]+shift[0])/shrink, (p[1]+shift[1])/shrink, (p[2]+shift[2])/shrink);
		}
	}
	return poly;
}

int main(int argc, char **argv){
	float shift[3] = {0,0,0};
	float shrink = 3;
	Polyhedron *poly = hispeed::read_polyhedron();
	Polyhedron apoly = adjust_polyhedron(shift,shrink,*poly);
	hispeed::write_polyhedron(&apoly, "/gisdata/shrinked.off");
	delete poly;
}
