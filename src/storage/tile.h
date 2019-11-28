/*
 * tile.h
 *
 *  Created on: Nov 15, 2019
 *      Author: teng
 */

#ifndef HISPEED_TILE_H_
#define HISPEED_TILE_H_
#include <stdio.h>

#include "himesh.h"

using namespace std;

namespace hispeed{


/*
 * each voxel contains the minimum boundary box
 * of a set of edges or triangles. It is an extension of
 * AAB with additional elements
 * */
class Voxel{
public:
	// point which the segments close wiht
	float core[3];
	// boundary box of the voxel
	aab box;
	// the pointer points to the segment data in this voxel
	float *data = NULL;
	int size = 0;
	HiMesh *mesh;
};

class Tile{
	aab box;
	FILE *dt_fs = NULL;
	std::vector<HiMesh_Wrapper *> objects;
	bool load();
	// retrieve the mesh of the voxel group with ID id on demand
	void retrieve_mesh(int id);
	void print();
public:

	Tile(std::string path);
	~Tile();

	HiMesh *get_mesh(int id, int lod){
		assert(id>=0&&id<objects.size());
		assert(lod>=0&&lod<=100);
		if(objects[id]->mesh==NULL){
			retrieve_mesh(id);
		}
		assert(objects[id]->mesh);
		objects[id]->mesh->advance_to(lod);
		return objects[id]->mesh;
	}
	HiMesh_Wrapper *get_mesh_wrapper(int id){
		assert(id>=0&&id<objects.size());
		return objects[id];
	}
	aab get_mbb(int id){
		assert(id>=0&&id<objects.size());
		return objects[id]->box;
	}
	int num_objects(){
		return objects.size();
	}
};

Tile *generate_tile();



}



#endif /* HISPEED_TILE_H_ */
