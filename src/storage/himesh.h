/*
 * himesh.h
 *
 *  Created on: Nov 28, 2019
 *      Author: teng
 */

#ifndef HIMESH_H_
#define HIMESH_H_

#include "../geometry/aab.h"
#include "../PPMC/mymesh.h"
#include "../PPMC/configuration.h"
#include "../spatial/spatial.h"

namespace hispeed{

class Voxel;
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

	// the buffer for filling segments to voxels
	float *segment_buffer = NULL;
public:
	HiMesh(const char *, long length);
	~HiMesh(){
		if(segment_buffer){
			delete segment_buffer;
		}
	}
	// added for HISPEED
	Skeleton *extract_skeleton();
	vector<Point> get_skeleton_points();
	vector<Voxel *> generate_voxels(int voxel_size);
	size_t fill_segments(float *segments);
	void fill_voxel(vector<Voxel *> &voxels);

	inline void get_vertices(std::vector<Point> &points){
		for(MyMesh::Vertex_iterator v = vertices_begin();
				v != vertices_end(); ++v){
			points.push_back(v->point());
		}
	}

	size_t size_of_edges(){
		return size_of_halfedges()/2;
	}
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
	vector<Voxel *> voxels;
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
	void writeMeshOff(){
		assert(mesh);
		stringstream ss;
		ss<<id<<".off";
		mesh->writeMeshOff(ss.str().c_str());
	}
	// fill the segments into voxels
	void fill_voxels(int lod){
		assert(mesh);
		mesh->advance_to(lod);
		mesh->fill_voxel(voxels);
	}

};

}




#endif /* HIMESH_H_ */
