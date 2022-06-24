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

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

// includes for the mesh simplification
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
// Visitor base
#include <CGAL/Surface_mesh_simplification/Edge_collapse_visitor_base.h>
// Extended polyhedron items which include an id() field
#include <CGAL/Polyhedron_items_with_id_3.h>
// Stop-condition policy
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>
// Non-default cost and placement policies
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>

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

namespace SMS = CGAL::Surface_mesh_simplification ;


namespace hispeed{

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
	size_t true_triangle_size();
	bool own_data = true;
	SegTree *aabb_tree = NULL;
	list<Segment> segments;
	std::map<int, Skeleton *> skeletons;
public:
	HiMesh(char *data, long length, bool own_data);
	HiMesh(char *data, long length);
	~HiMesh(){
		//release_buffer();
		this->clear_aabb_tree();
		for(std::map<int, Skeleton *>::iterator it=skeletons.begin();it!=skeletons.end();it++){
			delete it->second;
		}
		skeletons.clear();
	}

	aab get_box();

	Polyhedron *to_polyhedron();
	Skeleton *extract_skeleton();
	vector<Point> get_skeleton_points(int num_skeleton_points);
	vector<Voxel *> generate_voxels(int voxel_size);
	vector<Voxel *> voxelization(int voxel_size);
	string to_wkt();
	//float get_volume();

	size_t fill_segments(float *segments);
	size_t fill_triangles(float *triangles);
	size_t fill_voxels(vector<Voxel *> &voxels, enum data_type seg_or_triangle);
	size_t fill_vertices(float *&vertices);
	size_t fill_topology(unsigned short *&topology);
	void get_segments();
	SegTree *get_aabb_tree();
	void clear_aabb_tree();

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
void cgal_simplification(Polyhedron *poly, float ratio);


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
		//pthread_mutex_lock(&lock);
		mesh->advance_to(lod);
		//pthread_mutex_unlock(&lock);
	}
	// fill the segments into voxels
	// seg_tri: 0 for segments, 1 for triangle
	size_t fill_voxels(enum data_type seg_tri);
	size_t num_vertices(){
		return mesh->size_of_vertices();
	}

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
