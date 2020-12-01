/*
 * compress_single.cpp
 *
 *  Created on: Nov 30, 2020
 *      Author: teng
 */


#include "../util/util.h"
#include "../spatial/spatial.h"

using namespace std;

int main(int argc, char **argv){
	//MyMesh *mesh = hispeed::read_off(argv[1]);
	Polyhedron *poly = read_off_polyhedron(argv[1]);
}
