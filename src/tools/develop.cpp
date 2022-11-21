/*
 * develop.cpp
 *
 *  Created on: Nov 17, 2022
 *      Author: teng
 */

#include "../mymesh/mymesh.h"


using namespace hispeed;

int main(int argc, char **argv){


	Polyhedron *poly = new Polyhedron(0);
	poly->load(argv[1]);
	poly->remove_orphan_vertices();

	char path[256];
	sprintf(path,"/gisdata/mymesh/original.off");
	poly->dumpto(path);
	poly->evaluate();

	for(int i=0;i<10;i++){
		poly->reset_states();
		poly->compress();
		poly->print();
		sprintf(path,"/gisdata/mymesh/compressed_%d.off", i);
		poly->dumpto(path);
		poly->evaluate();
	}
	delete poly;

	return 0;
}





