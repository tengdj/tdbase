/*
 * tile.cpp
 *
 *  Created on: Nov 14, 2019
 *      Author: teng
 *
 *  the implementation of functions manipulating
 *  disk files and memory space which stores the
 *  compressed polyhedra and light weight index
 *  of each polyhedra. Each file contains the
 *  information for one single tile
 *
 */


#include "tile.h"

namespace hispeed{

// load meta data from file and construct the hierarchy structure
Tile::Tile(std::string path, size_t capacity){
	tile_path = path;
	tile_capacity = capacity;
	init();
}

Tile::~Tile(){
	for(HiMesh_Wrapper *h:objects){
		delete h;
	}
	if(data_buffer!=NULL){
		delete []data_buffer;
	}
}

// do the initialization job
void Tile::init(){
	struct timeval start = get_cur_time();
	if(!file_exist(tile_path.c_str())){
		log("%s does not exist", tile_path.c_str());
		exit(-1);
	}
	// load the meta data from the raw file or the cached meta file
	string meta_path = tile_path;
	boost::replace_all(meta_path, ".dt", ".mt");
	process_lock();
	// load the raw data into the buffer
	data_size = file_size(tile_path.c_str());
	data_buffer = new char[data_size];
	FILE *dt_fs = fopen(tile_path.c_str(), "r");
	assert(fread((void *)data_buffer, sizeof(char), data_size, dt_fs) == data_size);
	fclose(dt_fs);
	process_unlock();

	// parsing the metadata from the dt file
	size_t offset = 0;
	size_t index = 0;
	while(offset < data_size){
		size_t dsize = *(size_t *)(data_buffer + offset);
		offset += sizeof(size_t);

		// create a wrapper with the meta information
		HiMesh_Wrapper *w = new HiMesh_Wrapper(data_buffer + offset, dsize, index++);

		// next dsize bytes are for the mesh, skip them
		offset += dsize;
		// read the voxels into the wrapper
		size_t vnum = *(size_t *)(data_buffer + offset);
		offset += sizeof(size_t);

		for(int i=0;i<vnum;i++){
			Voxel *v = new Voxel();
			memcpy(v->low, data_buffer+offset, 3*sizeof(float));
			offset += 3*sizeof(float);
			memcpy(v->high, data_buffer+offset, 3*sizeof(float));
			offset += 3*sizeof(float);
			memcpy(v->core, data_buffer+offset, 3*sizeof(float));
			offset += 3*sizeof(float);

			w->voxels.push_back(v);
			w->box.update(*v);
		}
		objects.push_back(w);
		space.update(w->box);
	}

	// disable the feature of inner partitioning if configured
	if(!global_ctx.use_multimbb){
		for(auto *hmesh:objects){
			hmesh->disable_innerpart();
		}
	}
	logt("loaded %ld polyhedra in tile %s", start, objects.size(), tile_path.c_str());
}

HiMesh *Tile::get_mesh(int id){
	return objects[id]->mesh;
}

void Tile::decode_all(int lod){
	for(HiMesh_Wrapper *w:objects){
		w->decode_to(lod);
	}
}

OctreeNode *Tile::build_octree(size_t leaf_size){
	OctreeNode *octree = new OctreeNode(space, 0, leaf_size);
	for(HiMesh_Wrapper *w:objects){
		octree->addObject(&w->box);
	}
	return octree;
}

}

