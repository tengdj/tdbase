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
	himesh->advance_to(100);
	himesh->to_wkt();

	delete mesh;
	delete himesh;
}


