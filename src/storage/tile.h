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
#include "../PPMC/mymesh.h"
#include "../PPMC/configuration.h"
#include "../spatial/spatial.h"

using namespace std;

namespace hispeed{

class Tile;
class Voxel_group;

class Voxel{
public:
	int id = 0;
	// boundary box of the voxel
	aab box;
	// the decoded edges, or triangles in different LODs
	// the first value of each LOD represent the size of
	// the data for this LOD
	float *data[10];
	int size[10];

	// the voxel group current voxel belongs to
	Voxel_group *group = NULL;
	// can only be called by the Voxel_group class for cache evicting
	// the box and id will always exist for index looking up
	// but the detailed edges and triangles can be released if
	// memory is not big enough
	void clear();
	Voxel();
	~Voxel();

	int get_size(int lod);
	int get_total_size();
	void set_data(int lod, float *target_data);
	float *get_data(int lod);
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
public:

	int id = -1;
	MyMesh *mesh = NULL;
	char *mesh_data = NULL;
	// offset of the mesh in .dt file
	long offset = 0;
	// size of the mesh in .dt file
	long data_size = 0;

	// buffer for the decoded data shared by the voxels
	float *decoded_data[10];
	std::vector<Voxel *> voxels;
	// the tile containing this voxel group
	Tile *tile = NULL;


	Voxel_group();
	~Voxel_group();
	void add_voxel(aab box);

	// release the data for all voxels
	// this can be called when being evicted from cache
	void release_data();

	// decode the mesh to level of lod
	// also load the data from disk when needed
	void decode(int lod);

	bool persist(FILE *fs);
	bool load(FILE *fs);

};

class Tile{
	bool active = false;
	int id = -1;
	aab space;
	std::string meta_path;
	std::string data_path;
	std::vector<Voxel_group *> voxel_groups;
	FILE *dt_fs = NULL;

	// load the space and AAB of objects in this tile
	// can only be called by constructor
	bool load();

	void add_group(Voxel_group *group){
		group->id = this->voxel_groups.size();
		voxel_groups.push_back(group);
		group->tile = this;
	}
public:

	Tile(){}

	// load meta data from file
	// and construct the hierarchy structure
	// tile->voxel_groups->voxels
	Tile(std::string path);
	~Tile();
	void set_path(string path);
	bool is_active(){
		return active;
	}
	void add_polyhedron(MyMesh *mesh, long offset);
	bool persist();

	// retrieve the mesh of the voxel group with ID id on demand
	void retrieve_mesh(int id);

	// release the space for mesh object and the triangles
	// decoded in each voxel groups
	bool release_mesh(int id);

};

Tile *generate_tile();



}



#endif /* HISPEED_TILE_H_ */
