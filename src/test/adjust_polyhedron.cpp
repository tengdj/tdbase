/*
 * getoff.cpp
 *
 *  Created on: Nov 20, 2019
 *      Author: teng
 */

#include "../spatial/spatial.h"
#include "../spatial/himesh.h"

using namespace hispeed;

Polyhedron adjust_polyhedron(int shift[3], float shrink, Polyhedron &poly_o){

	Polyhedron poly;
	stringstream ss;
	ss << poly_o;
	ss >> poly;

	double min_co[3] = {DBL_MAX,DBL_MAX,DBL_MAX};
	double max_co[3] = {0,0,0};
	int counter = 0;
	int mean_co[3] = {0,0,0};
	for(Polyhedron::Vertex_iterator vi=poly.vertices_begin();vi!=poly.vertices_end();vi++){
		Point p = vi->point();
		for(int i=0;i<3;i++){
			if(min_co[i]>p[i]){
				min_co[i] = p[i];
			}
			if(max_co[i]<p[i]){
				max_co[i] = p[i];
			}
			mean_co[i] += p[i];
		}
		counter++;
	}

	for(int i=0;i<3;i++){
		mean_co[i] /= counter;
	}

	for(Polyhedron::Vertex_iterator vi=poly.vertices_begin();vi!=poly.vertices_end();vi++){
		Point p = vi->point();
		if(shrink>1){
			vi->point() = Point((p[0]+shift[0]-min_co[0])/shrink+mean_co[0],
							    (p[1]+shift[1]-min_co[1])/shrink+mean_co[1],
								(p[2]+shift[2]-min_co[2])/shrink+mean_co[2]);
		}else{
			vi->point() = Point((p[0]+shift[0])/shrink, (p[1]+shift[1])/shrink, (p[2]+shift[2])/shrink);
		}
	}
	return poly;
}

int main(int argc, char **argv){
	if(argc<=3){
		cout<<"usage: adjust shift shrink output_path"<<endl;
		return 0;
	}
	int sft = atoi(argv[1]);
	int shift[3] = {sft,sft,sft};
	float shrink = atoi(argv[2]);;
	Polyhedron *poly = hispeed::read_polyhedron();
	Polyhedron apoly = adjust_polyhedron(shift,shrink,*poly);
	hispeed::write_polyhedron(&apoly, argv[3]);
	delete poly;
}
