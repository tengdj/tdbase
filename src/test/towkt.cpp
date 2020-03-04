/*
 * towkt.cpp
 *
 *  Created on: Feb 15, 2020
 *      Author: teng
 */

#include "../spatial/spatial.h"
#include "../spatial/himesh.h"

using namespace hispeed;

int main(int argc, char **argv){

	MyMesh *mesh = hispeed::read_mesh();
	mesh->completeOperation();
	HiMesh *himesh = new HiMesh(mesh->p_data, mesh->dataOffset);
	int lod = 100;
	if(argc>2){
		lod = atoi(argv[2]);
	}
	himesh->advance_to(lod);
	himesh->to_wkt();

	delete mesh;
	delete himesh;
}


