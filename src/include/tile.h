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

namespace tdbase{

class Tile{
	aab space;
	std::vector<HiMesh_Wrapper *> objects;
	char *data_buffer = NULL;
	size_t data_size = 0;
	size_t tile_capacity = INT_MAX;
	string tile_path;

	OctreeNode *tree = NULL;
public:
	// for building tile instead of load from file
	Tile(std::vector<HiMesh_Wrapper *> &objs);
	Tile(std::string path, size_t capacity=LONG_MAX, bool active_load=true);
	~Tile();
	void load();

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

	HiMesh *get_mesh(int id);
	void decode_all(int lod = 100);
	OctreeNode *build_octree(size_t num_tiles);
	inline OctreeNode *get_octree(){
		return tree;
	}

	void dump_compressed(const char *path);
	void dump_raw(const char *path);

	void dump_sql(const char *path, const char *table);

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
