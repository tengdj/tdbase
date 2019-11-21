/*
 * getoff.cpp
 *
 *  Created on: Nov 20, 2019
 *      Author: teng
 */

#include "../spatial/spatial.h"

using namespace hispeed;

int main(int argc, char **argv){

	if(argc==1){
		cout<<"usage: getoff path/to/output"<<endl;
	}
	MyMesh *mesh = hispeed::read_mesh();
	mesh->writeMeshOff(argv[1]);

}
