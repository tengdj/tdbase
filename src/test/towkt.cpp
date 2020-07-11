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

	for(int i=0;i<1000;i++){
		MyMesh *mesh = hispeed::read_mesh();
		if(!mesh){
			break;
		}
		mesh->completeOperation();
		HiMesh *himesh = new HiMesh(mesh->p_data, mesh->dataOffset);
		int lod = 100;
		if(argc>2){
			lod = atoi(argv[1]);
		}
		himesh->advance_to(lod);
		himesh->to_wkt();

		delete mesh;
		delete himesh;
	}
}


