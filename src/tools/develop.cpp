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

	Vertex *v = poly->get_vertex(0);
	poly->remove_vertex(v);
	delete v;
	poly->print();
	poly->evaluate();
	delete poly;
	return 0;
}





