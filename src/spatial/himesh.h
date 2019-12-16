/*
 * himesh.h
 *
 *  Created on: Nov 28, 2019
 *      Author: teng
 */

#ifndef HIMESH_H_
#define HIMESH_H_

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_segment_primitive.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include "../geometry/aab.h"
#include "../geometry/geometry.h"
#include "../PPMC/mymesh.h"
#include "../PPMC/configuration.h"
#include "../spatial/spatial.h"



//CGAL::Polyhedron_3< MyKernel, MyItems > Polyhedron;

typedef MyKernel::FT FT;
typedef MyKernel::Segment_3 Segment;

typedef std::list<Segment>::iterator SegIterator;
typedef CGAL::AABB_segment_primitive<MyKernel, SegIterator> SegPrimitive;
typedef CGAL::AABB_traits<MyKernel, SegPrimitive> SegTraits;
typedef CGAL::AABB_tree<SegTraits> SegTree;

typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> TrianglePrimitive;
typedef CGAL::AABB_traits<MyKernel, TrianglePrimitive> TriangleTraits;
typedef CGAL::AABB_tree<TriangleTraits> TriangleTree;

namespace hispeed{

/*
 * each voxel contains the minimum boundary box
 * of a set of edges or triangles. It is an extension of
 * AAB with additional elements
 * */
class Voxel{
public:
	~Voxel(){
		reset();
	}
	// point which the segments close wiht
	float core[3];
	// boundary box of the voxel
	aab box;
	// the pointer and size of the segment/triangle data in this voxel
	map<int, float *> data;
	map<int, int> size;
	void reset(){
		for(map<int, float *>::iterator it=data.begin();it!=data.end();it++){
			if(it->second!=NULL){
				delete it->second;
				it->second = NULL;
			}
		}
		data.clear();
		size.clear();
	}
};

enum data_type{
	DT_Segment = 0,
	DT_Triangle
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
class HiMesh:public MyMesh{
	// the buffer for filling segments/triangles to voxels
	map<int, float *> segment_buffer;
	map<int, float *> triangle_buffer;
	size_t fill_segments(float *segments);
	size_t fill_triangles(float *triangles);
public:
	HiMesh(const char *, long length);
	~HiMesh(){
		release_buffer();
	}
	Polyhedron *to_polyhedron();
	void release_buffer();
	Skeleton *extract_skeleton();
	vector<Point> get_skeleton_points();
	vector<Voxel *> generate_voxels(int voxel_size);

	void fill_voxel(vector<Voxel *> &voxels, enum data_type seg_or_triangle);
	list<Segment> get_segments();
	SegTree *get_aabb_tree();

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
	bool filled = false;
	int id = -1;
	HiMesh *mesh = NULL;
	weighted_aab box;
	// used for retrieving compressed data from disk
	size_t offset = 0;
	size_t data_size = 0;
	pthread_mutex_t lock;
	HiMesh_Wrapper(){
		pthread_mutex_init(&lock, NULL);
	}
	~HiMesh_Wrapper(){
		for(Voxel *v:voxels){
			delete v;
		}
		voxels.clear();
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
	void advance_to(int lod){
		assert(mesh);
		pthread_mutex_lock(&lock);
		mesh->advance_to(lod);
		pthread_mutex_unlock(&lock);
	}
	// fill the segments into voxels
	// seg_tri: 0 for segments, 1 for triangle
	void fill_voxels(enum data_type seg_tri);

	void reset(){
		pthread_mutex_lock(&lock);
		for(Voxel *v:voxels){
			v->reset();
		}
		pthread_mutex_unlock(&lock);
	}


};

TriangleTree *get_aabb_tree(Polyhedron *p);


}

#endif /* HIMESH_H_ */
