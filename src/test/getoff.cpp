/*
 * getoff.cpp
 *
 *  Created on: Nov 20, 2019
 *      Author: teng
 */

#include "../spatial/spatial.h"

using namespace hispeed;

int main(int argc, char **argv){

	if(argc<2){
		log("usage: getoff path/to/output");
		return 0;
	}
	MyMesh *mesh = hispeed::read_mesh();
	mesh->writeMeshOff(argv[1]);
	delete mesh;
}
