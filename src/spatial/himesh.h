/*
 * himesh.h
 *
 *  Created on: Nov 28, 2019
 *      Author: teng
 */

#ifndef HIMESH_H_
#define HIMESH_H_

#define CGAL_EIGEN3_ENABLED


#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/box_intersection_d.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/intersections.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/convex_decomposition_3.h>
#include <CGAL/Tetrahedron_3.h>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/extract_mean_curvature_flow_skeleton.h>
#include <CGAL/boost/graph/split_graph_into_polylines.h>
#include <boost/foreach.hpp>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_segment_primitive.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

// includes for the mesh simplification
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Edge_collapse_visitor_base.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>

#include <boost/algorithm/string/replace.hpp>

#include "../geometry/aab.h"
#include "../geometry/geometry.h"
#include "../PPMC/mymesh.h"
#include "../PPMC/configuration.h"

// templates for AABB tree
typedef MyKernel::FT FT;
typedef MyKernel::Segment_3 Segment;
typedef MyKernel::Triangle_3 Triangle;
typedef std::list<Segment>::iterator SegIterator;
typedef CGAL::AABB_segment_primitive<MyKernel, SegIterator> SegPrimitive;
typedef CGAL::AABB_traits<MyKernel, SegPrimitive> SegTraits;
typedef CGAL::AABB_tree<SegTraits> SegTree;

typedef std::list<Triangle>::iterator Iterator;
typedef CGAL::AABB_triangle_primitive<MyKernel, Iterator> TrianglePrimitive;
typedef CGAL::AABB_traits<MyKernel, TrianglePrimitive> TriangleTraits;
typedef CGAL::AABB_tree<TriangleTraits> TriangleTree;

//
typedef CGAL::Polyhedron_3<MyKernel, CGAL::Polyhedron_items_with_id_3> Polyhedron;
typedef boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;

typedef CGAL::Surface_mesh<Point>                             Triangle_mesh;
typedef CGAL::Mean_curvature_flow_skeletonization<Triangle_mesh> Skeletonization;
typedef Skeletonization::Skeleton                             Skeleton;
typedef Skeleton::vertex_descriptor                           Skeleton_vertex;
typedef Skeleton::edge_descriptor                             Skeleton_edge;

typedef CGAL::Delaunay_triangulation_3<MyKernel, CGAL::Fast_location> Delaunay;

typedef CGAL::Triangulation_3<MyKernel> Triangulation;
typedef Triangulation::Tetrahedron 	Tetrahedron;

namespace SMS = CGAL::Surface_mesh_simplification ;
//
//typedef CGAL::Nef_polyhedron_3<MyKernel> Nef_polyhedron;
//typedef Nef_polyhedron::Volume_const_iterator Volume_const_iterator;
//typedef Nef_polyhedron::Shell_entry_const_iterator Shell_entry_const_iterator;

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
	SegTree *segment_tree = NULL;
	TriangleTree *triangle_tree = NULL;
	list<Segment> segments;
	list<Triangle> triangles;
	void get_segments();
	void get_triangles();
public:
	HiMesh(char *data, long length, bool own_data);
	HiMesh(char *data, long length);
	~HiMesh();

	aab get_box();

	Polyhedron *to_polyhedron();
	Polyhedron *to_triangulated_polyhedron();
	Skeleton *extract_skeleton();
	vector<Point> get_skeleton_points(int num_skeleton_points);
	vector<Voxel *> generate_voxels_skeleton(int voxel_size);
	vector<Voxel *> voxelization(int voxel_size);
	string to_wkt();
	float get_volume();

	size_t fill_segments(float *segments);
	size_t fill_triangles(float *triangles);
	size_t fill_voxels(vector<Voxel *> &voxels, enum data_type seg_or_triangle);
	size_t fill_vertices(float *&vertices);
	size_t fill_topology(unsigned short *&topology);
	SegTree *get_aabb_tree_segment();
	TriangleTree *get_aabb_tree_triangle();

	void clear_aabb_tree();

	void get_vertices(std::vector<Point> &points);
	size_t size_of_edges();
	void advance_to(int lod);
};

/*
 * a wrapper for describing a mesh. for some cases we may not need to truly
 * parse the mesh out from disk, but use the bounding box is good enough.
 * */
class HiMesh_Wrapper{
public:
	vector<Voxel *> voxels;
	bool filled = false;
	size_t id = -1;
	HiMesh *mesh = NULL;
	weighted_aab box;
	// used for retrieving compressed data from disk
	size_t offset = 0;
	size_t data_size = 0;

	pthread_mutex_t lock;

	int candidate_confirmed = 0;
public:
	HiMesh_Wrapper();
	~HiMesh_Wrapper();
	void writeMeshOff();
	void advance_to(int lod);
	// fill the segments into voxels
	// seg_tri: 0 for segments, 1 for triangle
	size_t fill_voxels(enum data_type seg_tri);

	size_t num_vertices();
	void reset();

};

// some general utility functions

void cgal_simplification(Polyhedron *poly, float ratio);

Polyhedron *make_cube(aab box);
Polyhedron *make_cubes(vector<aab *> &boxes);

void write_box(aab box, int id, string prefix="");
void write_box(aab box, const char *path);
void write_polyhedron(Polyhedron *mesh, const char *path);
void write_polyhedron(Polyhedron *mesh, int id);
void write_voxels(vector<Voxel *> voxels, const char *path);
void write_points(vector<Point> &skeleton, const char *path);
string read_off_stdin();
string polyhedron_to_wkt(Polyhedron *poly);

// some utility functions to operate mesh polyhedrons
extern MyMesh *get_mesh(string input, bool complete_compress = false);
extern MyMesh *read_mesh();
extern MyMesh *read_off(char *path);
extern Polyhedron *read_off_polyhedron(char *path);

extern MyMesh *decompress_mesh(MyMesh *compressed, int lod, bool complete_operation = false);
extern MyMesh *decompress_mesh(char *data, size_t length, bool complete_operation = false);

float get_volume(Polyhedron *polyhedron);
Polyhedron *read_polyhedron();
Polyhedron *read_polyhedron(string &str);
Polyhedron adjust_polyhedron(int shift[3], float shrink, Polyhedron *poly_o);

inline float distance(Point p1, Point p2){
	float dist = 0;
	for(int i=0;i<3;i++){
		dist += (p2[i]-p1[i])*(p2[i]-p1[i]);
	}
	return dist;
}

// manhattan distance
inline float mdistance(Point p1, Point p2){
	float dist = 0;
	for(int i=0;i<3;i++){
		dist += abs((float)(p2[i]-p1[i]));
	}
	return dist;
}


}

#endif /* HIMESH_H_ */
