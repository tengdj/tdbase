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

#include "aab.h"
#include "geometry.h"
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
	SegTree *segment_tree = NULL;
	TriangleTree *triangle_tree = NULL;
	list<Segment> segments;
	list<Triangle> triangles;
	list<Point> vertices;
public:
	HiMesh(char *data, size_t dsize);
	HiMesh(MyMesh *mesh);
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

	// get the elements
	size_t fill_segments(float *&segments);
	size_t fill_triangles(float *&triangles);
	size_t fill_hausdorf_distances(float *&hausdorf);
	pair<float, float> get_triangle_hausdorf(int tri_id = -1);
	size_t fill_voxels(vector<Voxel *> &voxels);
	size_t fill_vertices(float *&vertices);
	list<Segment> get_segments();
	list<Triangle> get_triangles();
	list<Point> get_vertices();

	// query
	SegTree *get_aabb_tree_segment();
	TriangleTree *get_aabb_tree_triangle();
	float distance(HiMesh *target);
	float distance_tree(HiMesh *target);
	bool intersect(HiMesh *target);
	bool intersect_tree(HiMesh *target);

	void clear_aabb_tree();

	size_t size_of_edges();

	void advance_to(int lod);

	// validation
	bool has_same_vertices();
};

/*
 * a wrapper for describing a mesh. for some cases we may not need to truly
 * parse the mesh out from disk, but use the bounding box is good enough.
 * */
class HiMesh_Wrapper{
public:
	// description of the polyhedron
	size_t id = -1;
	weighted_aab box;
	vector<Voxel *> voxels;

	// used for retrieving compressed data from disk
	HiMesh *mesh = NULL;
	size_t offset = 0;
	size_t data_size = 0;

	pthread_mutex_t lock;
	vector<HiMesh_Wrapper *> results;
public:
	HiMesh_Wrapper();
	~HiMesh_Wrapper();
	void writeMeshOff();
	void advance_to(int lod);
	// fill the triangles into voxels
	size_t fill_voxels();

	size_t num_vertices();
	void reset();
	void report_result(HiMesh_Wrapper *result);

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
extern MyMesh *parse_mesh(string input, bool complete_compress = false);
extern MyMesh *read_mesh();
extern MyMesh *read_mesh(char *path);

extern MyMesh *decompress_mesh(MyMesh *compressed, int lod, bool complete_operation = false);

float get_volume(Polyhedron *polyhedron);
Polyhedron *read_polyhedron();
extern Polyhedron *read_polyhedron(const char *path);
vector<Polyhedron *> read_polyhedrons(const char *path, size_t load_num = LONG_MAX);

Polyhedron *parse_polyhedron(string &str);
Polyhedron adjust_polyhedron(int shift[3], float shrink, Polyhedron *poly_o);

}

#endif /* HIMESH_H_ */
