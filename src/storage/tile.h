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
#include <pthread.h>

using namespace std;

namespace hispeed{

class Tile{
	pthread_mutex_t read_lock;
	aab box;
	std::vector<HiMesh_Wrapper *> objects;
	FILE *dt_fs = NULL;
	bool load(string path, int max_objects=INT_MAX);
	bool persist(string path);
	bool parse_raw();
	// retrieve the data of the mesh with ID id on demand
	void retrieve_mesh(int id);

public:
	// for building tile instead of load from file
	Tile(){};
	void add_raw(char *data);
	Tile(std::string path, size_t capacity=LONG_MAX);
	~Tile();
	void disable_innerpart();

	// for profiling performance
	double decode_time = 0;
	double retrieve_time = 0;
	double advance_time = 0;
	double disk_time = 0;
	double malloc_time = 0;
	double newmesh_time = 0;
	void reset_time(){
		decode_time = 0;
		retrieve_time = 0;
		advance_time = 0;
		disk_time = 0;
		malloc_time = 0;
		newmesh_time = 0;
	}

	void print_time(){
		cerr<<"\ndecoding time\t"<<decode_time
			<<"\n\tretrieve time\t"<< retrieve_time
			<<"\n\t\tdisk time\t" << disk_time
			<<"\n\t\tmalloc time\t"<< malloc_time
			<<"\n\t\tnewmesh time\t"<< newmesh_time
			<<"\n\tadvance time\t"<< advance_time
			<<endl<<endl;
	}

	void decode_to(int id, int lod);
	HiMesh_Wrapper *get_mesh_wrapper(int id){
		assert(id>=0&&id<objects.size());
		return objects[id];
	}
	aab get_mbb(int id){
		assert(id>=0&&id<objects.size());
		return objects[id]->box;
	}
	HiMesh *get_mesh(int id){
		assert(get_mesh_wrapper(id)->mesh && "the mesh must be retrieved before can be returned");
		return get_mesh_wrapper(id)->mesh;
	}
	size_t num_objects(){
		return objects.size();
	}

	void retrieve_all(){
		for(HiMesh_Wrapper *w:objects){
			retrieve_mesh(w->id);
		}
	}

	void advance_all(int lod){
		retrieve_all();
		for(HiMesh_Wrapper *w:objects){
			w->advance_to(lod);
		}
	}

	OctreeNode *build_octree(size_t num_tiles);
	//SpatialIndex::ISpatialIndex *build_rtree();

};

}



#endif /* HISPEED_TILE_H_ */
