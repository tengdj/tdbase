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
#include "index.h"
#include <pthread.h>

using namespace std;

namespace hispeed{

class Tile{
	aab box;
	std::vector<HiMesh_Wrapper *> objects;
	char *data_buffer = NULL;
	size_t tile_capacity = INT_MAX;
	string tile_path;

	void init();
	bool load(string meta_path, int max_objects=INT_MAX);
	bool persist(string meta_path);
	bool parse_raw();
	// retrieve the data of the mesh with ID id on demand
	void retrieve_mesh(size_t id);
public:
	// for building tile instead of load from file
	Tile(std::string path, size_t capacity=LONG_MAX);
	~Tile();

	inline HiMesh_Wrapper *get_mesh_wrapper(int id){
		assert(id>=0&&id<objects.size());
		return objects[id];
	}
	inline aab get_mbb(int id){
		assert(id>=0&&id<objects.size());
		return objects[id]->box;
	}

	inline size_t num_objects(){
		return objects.size();
	}

	void decode_to(size_t id, uint lod);
	HiMesh *get_mesh(int id);
	void retrieve_all();
	void decode_all(int lod = 100);
	char *retrieve_data(int id);
	size_t get_object_data_size(int id);
	OctreeNode *build_octree(size_t num_tiles);

	// for profiling performance
private:
	double decode_time = 0;
	double retrieve_time = 0;
	double advance_time = 0;
	double newmesh_time = 0;
public:
	void reset_time(){
		decode_time = 0;
		retrieve_time = 0;
		advance_time = 0;
		newmesh_time = 0;
	}
	void print_time(){
		cerr<<"\ndecoding time\t"<<decode_time
			<<"\n\tretrieve time\t"<< retrieve_time
			<<"\n\t\tnewmesh time\t"<< newmesh_time
			<<"\n\tadvance time\t"<< advance_time
			<<endl<<endl;
	}
};

}



#endif /* HISPEED_TILE_H_ */
