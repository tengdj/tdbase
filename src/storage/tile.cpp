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

// retrieve the mesh of the voxel group with ID id on demand
void Tile::retrieve_mesh(size_t id){
	assert(id>=0 && id<objects.size());
	HiMesh_Wrapper *wrapper = objects[id];
	if(wrapper->mesh==NULL){
		timeval cur = hispeed::get_cur_time();
		wrapper->mesh = new HiMesh(data_buffer+wrapper->offset, wrapper->data_size);
		wrapper->mesh->id = id;
		newmesh_time += hispeed::get_time_elapsed(cur, true);
	}
	assert(wrapper->mesh);
}

HiMesh *Tile::get_mesh(int id){
	if(!get_mesh_wrapper(id)->mesh){
		retrieve_mesh(id);
	}
	assert(get_mesh_wrapper(id)->mesh && "the mesh must be retrieved before can be returned");
	return get_mesh_wrapper(id)->mesh;
}


void Tile::retrieve_all(){
	for(HiMesh_Wrapper *w:objects){
		retrieve_mesh(w->id);
	}
}

void Tile::decode_all(int lod){
	retrieve_all();
	for(HiMesh_Wrapper *w:objects){
		w->decode_to(lod);
	}
}

char *Tile::retrieve_data(int id){
	char *ret = new char[objects[id]->data_size];
	memcpy(ret, data_buffer+objects[id]->offset, objects[id]->data_size);
	return ret;
}

size_t Tile::get_object_data_size(int id){
	return objects[id]->data_size;
}

OctreeNode *Tile::build_octree(size_t leaf_size){
	OctreeNode *octree = new OctreeNode(box, 0, leaf_size);
	for(HiMesh_Wrapper *w:objects){
		octree->addObject(&w->box);
	}
	return octree;
}


void Tile::decode_to(size_t id, uint lod){
	assert(lod>=0 && lod<=100);
	assert(id>=0&&id<objects.size());
	timeval cur = hispeed::get_cur_time();
	timeval start = hispeed::get_cur_time();
	retrieve_mesh(id);
	retrieve_time += hispeed::get_time_elapsed(cur,true);
	objects[id]->decode_to(lod);
	advance_time += hispeed::get_time_elapsed(cur,true);
	decode_time += hispeed::get_time_elapsed(start,true);
}

}

