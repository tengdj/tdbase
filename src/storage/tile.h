/*
 * tile.h
 *
 *  Created on: Nov 15, 2019
 *      Author: teng
 */

#ifndef HISPEED_TILE_H_
#define HISPEED_TILE_H_
#include <stdio.h>

#include "../spatial/himesh.h"
#include "../index/index.h"

using namespace std;

namespace hispeed{

class Tile{
	aab box;
	std::vector<HiMesh_Wrapper *> objects;
	FILE *dt_fs = NULL;
	bool load();
	// retrieve the mesh of the voxel group with ID id on demand
	void retrieve_mesh(int id);
public:

	Tile(std::string path);
	~Tile();

	HiMesh *get_mesh(int id, int lod);
	HiMesh_Wrapper *get_mesh_wrapper(int id){
		assert(id>=0&&id<objects.size());
		return objects[id];
	}
	aab get_mbb(int id){
		assert(id>=0&&id<objects.size());
		return objects[id]->box.box;
	}
	int num_objects(){
		return objects.size();
	}

	OctreeNode *build_octree(int num_tiles);


};

}



#endif /* HISPEED_TILE_H_ */
