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
		size = 0;
		data = NULL;
	}
	// point which the segments close wiht
	float core[3];
	// boundary box of the voxel
	aab box;
	// the pointer points to the segment data in this voxel
	float *data = NULL;
	int size = 0;
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

	// the buffer for filling segments to voxels
	float *data_buffer = NULL;
public:
	HiMesh(const char *, long length);
	~HiMesh(){
		release_buffer();
	}
	void release_buffer(){
		if(data_buffer){
			delete data_buffer;
			data_buffer = NULL;
		}
	}
	// added for HISPEED
	Skeleton *extract_skeleton();
	vector<Point> get_skeleton_points();
	vector<Voxel *> generate_voxels(int voxel_size);
	size_t fill_segments(float *segments);
	size_t fill_triangles(float *triangles);

	void fill_voxel(vector<Voxel *> &voxels, int seg_or_triangle);
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
	aab box;
	// used for retrieving compressed data from disk
	long offset = 0;
	long data_size = 0;
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
	// fill the segments into voxels
	// seg_tri: 0 for segments, 1 for triangle
	void fill_voxels(int lod, int seg_tri){
		assert(mesh);
		mesh->advance_to(lod);
		mesh->fill_voxel(voxels, seg_tri);
	}

	void reset(){
		for(Voxel *v:voxels){
			v->data = NULL;
			v->size = 0;
		}
		if(mesh){
			mesh->release_buffer();
		}
	}

};

TriangleTree *get_aabb_tree(Polyhedron *p);



}

#endif /* HIMESH_H_ */
