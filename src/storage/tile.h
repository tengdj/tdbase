/*
 * tile.h
 *
 *  Created on: Nov 15, 2019
 *      Author: teng
 */

#ifndef HISPEED_TILE_H_
#define HISPEED_TILE_H_
#include <stdio.h>

#include "../geometry/aab.h"

namespace hispeed{

class Voxel_group;
class Voxel{

	long id = 0;
	// boundary box of the voxel
	aab box;
	// the decoded edges, or triangles in different LODs
	// the first value of each LOD represent the size of
	// the data for this LOD
	char *data[10];
	int size = 0;

	// the voxel group current voxel belongs to
	Voxel_group *group = NULL;
	friend class voxel_group;
protected:
	// can only be called by the MyMysh class for cache evicting
	// the box and id will always exist for index looking up
	// but the detailed edges and triangles can be released if
	// memory is not big enough
	void clear(){
		for(int i=0;i<10;i++){
			if(data[i]!=NULL){
				delete data[i];
				data[i] = NULL;
			}
		}
		size = 0;
	}
public:
	Voxel(){
		for(int i=0;i<10;i++){
			data[i] = NULL;
		}
	};
	~Voxel(){
		clear();
	}

	int get_size(){
		return size;
	}

	float *get_edges(int lod){
		assert(lod>=0&&lod<10);
		if(data[lod]==NULL){
			group->decode(lod);
		}
		return (float *)data[lod];
	}

	float *get_triangles(int lod){
		assert(lod>=0&&lod<10);
		return (float *)data[lod];
	}

};


/*
 *
 * each voxel group is mapped to an independent polyhedron
 * it contains a list of voxels, and initially being assigned
 * with the ABB of the voxel, for indexing. The true surfaces
 * and edges for different levels of details will be decoded and
 * released on-demand. The spaces taken by them is managed by
 * a global cache.
 *
 * */
class Voxel_group{

	// the file containing the compressed
	// meshes
	FILE* fs = NULL;
	MyMesh *mesh = NULL;
	// offset of the mesh in file
	long offset;
	// size of the mesh
	long data_size;
	std::vector<Voxel *> voxels;
	long num_vexels;

public:

	long get_data_size(){
		return data_size;
	}
	Voxel_group(FILE *fs, long offset, long data_length){
		this->fs = fs;
		this->offset = offset;
		this->data_size = size;
		num_vexels = 0;
	}

	~Voxel_group(){
		if(mesh!=NULL){
			delete mesh;
		}
		for(Voxel *v:voxels){
			if(v!=NULL){
				delete v;
			}
		}
		voxels.clear();
	}

	// release the data for all voxels
	// this can be called when being evicted from cache
	void release_data();

	// decode the mesh to level of lod
	// also load the data from disk when needed
	void decode(int lod);

};

class Tile{
	bool active = false;
	long id;
	aab space;
	std::string file_path;
	std::vector<Voxel_group *> voxels;
	FILE* fs = NULL;
public:

	Tile(){
		id = -1;
	}

	//load from file
	Tile(std::string path){
		if(load()){
			active = true;
		}
	}


	~Tile(){
		for(voxel *v:voxels){
			delete v;
		}
		if(fs!=NULL){
			flose(fs);
			fs = NULL;
		}
	}

	// load the space and AAB of objects in this tile
	bool load();
	// retrieve edges or triangles for certain voxel
	// to memory space given by data
	bool retrieve(int voxel_id, float *data);

	// persist the content of the tile to disk file
	bool persist();


};


}



#endif /* HISPEED_TILE_H_ */
