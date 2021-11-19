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
	stringstream ss;
	ss<<*poly;

	timeval start = hispeed::get_cur_time();
	MyMesh *mesh = hispeed::get_mesh(ss.str());
	logt("get mesh",start);
	mesh->completeOperation();
	logt("encoding takes",start);
	HiMesh *himesh = new HiMesh(mesh->p_data, mesh->dataOffset);
	himesh->advance_to(atoi(argv[1]));
	logt("decoding takes",start);
	cout<<himesh->size_of_vertices()<<" "<<himesh->get_skeleton_points(150).size()<<endl;
	logt("extract skeleton",start);
	vector<Voxel *> voxels = himesh->generate_voxels(30);
	logt("generate voxels",start);
	himesh->get_aabb_tree();
	logt("building index",start);
	for(int j=0;j<voxels.size();j++){
		hispeed::write_polyhedron(hispeed::make_cube(voxels[j]->box),j);
	}
	himesh->writeMeshOff("v.off");



//	for(int i=0;i<=100;i+=10){
//		himesh->advance_to(i);
//		cout<<i<<","<<himesh->size_of_vertices()<<","<<himesh->size_of_edges()<<endl;
//	}
//	cout<<endl;
//	for(int i=1;i<=1000;i+=1){
//		vector<Voxel *> boxes = himesh->generate_voxels(i);
//		int num_boxes = boxes.size();
//		float size = 0;
//		if(i==100){
//			for(int j=0;j<boxes.size();j++){
//				hispeed::write_polyhedron(hispeed::make_cube(boxes[j]->box),j);
//			}
//		}
//		for(Voxel *v:boxes){
//			size += v->box.volume();
//			delete v;
//		}
//		boxes.clear();
//		cout<<i<<","<<num_boxes<<","<<size/volume<<endl;
//
//	}



	delete himesh;
	delete mesh;

}
