/*
 * skeleton.cpp
 *
 *  Created on: Nov 12, 2019
 *      Author: teng
 */


#include "../spatial/spatial.h"
#include "../PPMC/ppmc.h"
#include "../util/util.h"

using namespace std;
using namespace hispeed;
using namespace CGAL;

int main(int argc, char **argv){

	MyMesh *mesh = read_mesh();
	mesh->set_skeleton_sample_rate(20);
	mesh->generate_mbbs();

	std::vector<aab> skeleton_mbbs = mesh->get_mbbs();

	int index = 0;
	char path[256];
	for(aab box:skeleton_mbbs){
		Polyhedron *pbox = make_cube(box);
		std::stringstream os;
		sprintf(path,"offs/%d.off",index++);
		write_polyhedron(pbox, path);
		delete pbox;
	}

	skeleton_mbbs.clear();
	delete mesh;
}



