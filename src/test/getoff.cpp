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
		cout<<"usage: getoff path/to/output"<<endl;
		return 0;
	}
	MyMesh *mesh = hispeed::read_mesh();
	mesh->writeMeshOff(argv[1]);
	cout<<*mesh<<endl;
	delete mesh;
}
