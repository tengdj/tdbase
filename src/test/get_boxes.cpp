/*
 * joiner.cpp
 *
 *  Created on: Nov 20, 2019
 *      Author: teng
 */

#include <boost/program_options.hpp>
#include <vector>

#include "../join/SpatialJoin.h"
#include "../spatial/himesh.h"
#include "../index/index.h"

using namespace std;
using namespace hispeed;
namespace po = boost::program_options;

int main(int argc, char **argv){
	struct timeval start = get_cur_time();
	MyMesh *mesh = hispeed::read_mesh();
	mesh->completeOperation();
	HiMesh *himesh = new HiMesh(mesh->p_data, mesh->dataOffset);
	himesh->advance_to(100);
	int sample_rate = 100;
	if(argc>=2){
		sample_rate = atoi(argv[1]);
	}
	himesh->writeMeshOff("/gisdata/vessel.off");

	vector<Voxel *> voxels = himesh->generate_voxels(sample_rate);

	aab b;
	for(int i=0;i<voxels.size();i++){
		b.update(voxels[i]->box);
		hispeed::write_box(voxels[i]->box, i, "/gisdata/boxes");
	}
	hispeed::write_box(b, "/gisdata/aab.off");
	delete himesh;
}
