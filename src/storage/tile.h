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
class HiMesh:public MyMesh{
	// the tile containing this voxel group
	std::vector<Voxel *> voxels;
	float *segment_buffer = NULL;
public:
	HiMesh(const char *, long length);
	~HiMesh(){
		if(segment_buffer){
			delete segment_buffer;
			voxels.clear();
		}
	}
	// added for HISPEED
	Skeleton *extract_skeleton();
	// generate local minimum boundary boxes
	// with the skeleton extracted
	vector<aab> generate_mbbs(int voxel_size);

	inline void get_vertices(std::vector<Point> &points){
		for(MyMesh::Vertex_iterator v = vertices_begin();
				v != vertices_end(); ++v){
			points.push_back(v->point());
		}
	}

	int get_segment_num(){
		return size_of_halfedges()/2;
	}
	void get_segments(float *segments=NULL);
	void advance_to(int lod);
};

/*
 * a wrapper for describing a mesh.
 * for some cases we may not need to truly
 * parse the mesh out from disk, but use the boundary box
 * is good enough.
 * */
class HiMesh_Wrapper{
public:
	int id = -1;
	HiMesh *mesh = NULL;
	aab box;
	// used for retrieving compressed data from disk
	long offset = 0;
	long data_size = 0;
	~HiMesh_Wrapper(){
		if(mesh){
			delete mesh;
		}
	}
};

class Tile{
	aab box;
	std::string prefix;
	FILE *dt_fs = NULL;

	std::vector<HiMesh_Wrapper *> objects;

	bool load();
	bool persist();
	bool parse_raw();
	// retrieve the mesh of the voxel group with ID id on demand
	void retrieve_mesh(int id);

	void print();
	std::string get_meta_path(){
		return prefix+".mt";
	}
	std::string get_data_path(){
		return prefix+".dt";
	}
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
