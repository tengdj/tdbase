/*
 * getoff.cpp
 *
 *  Created on: Nov 20, 2019
 *      Author: teng
 */

#include "../spatial/spatial.h"
#include "../spatial/himesh.h"

using namespace hispeed;

int main(int argc, char **argv){

//	if(argc<2){
//		log("usage: getoff path/to/output");
//		return 0;
//	}
	Polyhedron *poly = hispeed::read_polyhedron();
	float volume = -CGAL::Polygon_mesh_processing::volume(*poly);
	delete poly;
	MyMesh *mesh = hispeed::read_mesh();
	mesh->completeOperation();
	HiMesh *himesh = new HiMesh(mesh->p_data, mesh->dataOffset);
	for(int i=0;i<=100;i+=10){
		himesh->advance_to(i);
		cout<<i<<","<<himesh->size_of_vertices()<<","<<himesh->size_of_edges()<<endl;
	}
	for(int i=370;i<=1000;i+=5){
		vector<Voxel *> boxes = himesh->generate_voxels(i);
		int num_boxes = boxes.size();
		float size = 0;
		for(Voxel *v:boxes){
			size += v->box.volume();
			delete v;
		}
		boxes.clear();
		cout<<i<<","<<num_boxes<<","<<size/volume<<endl;
	}
	delete mesh;
	delete himesh;

}
