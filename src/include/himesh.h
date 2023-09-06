/*****************************************************************************
* Copyright (C) 2011 Adrien Maglo and Clément Courbet
*
* This file is part of PPMC.
*
* PPMC is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* PPMC is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with PPMC.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/

#ifndef PROGRESSIVEPOLYGONS_MYMESH_H
#define PROGRESSIVEPOLYGONS_MYMESH_H

#ifndef CGAL_EIGEN3_ENABLED
#define CGAL_EIGEN3_ENABLED
#endif

#include <iostream>
#include <fstream>
#include <stdint.h>
#include <queue>
#include <assert.h>
#include <utility>
#include <unordered_map>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/circulator.h>
#include <CGAL/bounding_box.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

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


#include <CGAL/config.h>
#include <boost/algorithm/string/replace.hpp>


#include "aab.h"
#include "util.h"
#include "geometry.h"
#include "query_context.h"

// Size of the compressed data buffer.
// for configuration
#define BUFFER_SIZE 10 * 1024 * 1024

#define NB_BITS_FACE_DEGREE_BASE_MESH 3

const int COMPRESSION_MODE_ID = 0;
const int DECOMPRESSION_MODE_ID = 1;

#define INV_ALPHA 2
#define INV_GAMMA 2

#define PPMC_RANDOM_CONSTANT 0315

using namespace std;
namespace SMS = CGAL::Surface_mesh_simplification ;

// definition for the CGAL library
//typedef CGAL::Exact_predicates_exact_constructions_kernel MyKernel;
typedef CGAL::Simple_cartesian<float> MyKernel;

typedef MyKernel::Point_3 Point;
typedef MyKernel::Vector_3 Vector;

typedef CGAL::Simple_cartesian<double> MyKernelDouble;
typedef MyKernelDouble::Vector_3 VectorDouble;

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

//typedef CGAL::Nef_polyhedron_3<MyKernel> Nef_polyhedron;
//typedef Nef_polyhedron::Volume_const_iterator Volume_const_iterator;
//typedef Nef_polyhedron::Shell_entry_const_iterator Shell_entry_const_iterator;

typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<MyKernel, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef Tree::Point_and_primitive_id Point_and_primitive_id;

namespace hispeed{

// the builder for reading the base mesh
template <class HDS> class MyMeshBaseBuilder : public CGAL::Modifier_base<HDS>
{
public:
    MyMeshBaseBuilder(std::deque<Point> *p_pointDeque, std::deque<uint32_t *> *p_faceDeque)
        : p_pointDeque(p_pointDeque), p_faceDeque(p_faceDeque) {}

    void operator()(HDS& hds)
    {
        CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);

        size_t nbVertices = p_pointDeque->size();
        size_t nbFaces = p_faceDeque->size();

        B.begin_surface(nbVertices, nbFaces);

        for (unsigned i = 0; i < nbVertices; ++i)
            B.add_vertex(p_pointDeque->at(i));

        for (unsigned i = 0; i < nbFaces; ++i)
        {
            B.begin_facet();
            uint32_t *f = p_faceDeque->at(i);
            for (unsigned j = 1; j < f[0] + 1; ++j)
                B.add_vertex_to_facet(f[j]);
            B.end_facet();
        }

        B.end_surface();
    }

private:
    std::deque<Point> *p_pointDeque;
    std::deque<uint32_t *> *p_faceDeque;
};

// My vertex type has a isConquered flag
template <class Refs>
class MyVertex : public CGAL::HalfedgeDS_vertex_base<Refs,CGAL::Tag_true, Point>
{
    enum Flag {Unconquered=0, Conquered=1};

	Flag flag = Unconquered;
	unsigned int id = 0;
  public:
    MyVertex(): CGAL::HalfedgeDS_vertex_base<Refs,CGAL::Tag_true, Point>(){}
	MyVertex(const Point &p): CGAL::HalfedgeDS_vertex_base<Refs,CGAL::Tag_true, Point>(p){}

	inline void resetState()
	{
	  flag=Unconquered;
	}

	inline bool isConquered() const
	{
	  return flag==Conquered;
	}

	inline void setConquered()
	{
	  flag=Conquered;
	}

	inline size_t getId() const
	{
		return id;
	}

	inline void setId(size_t nId)
	{
		id = nId;
	}
};

template <class Refs>
class MyHalfedge : public CGAL::HalfedgeDS_halfedge_base<Refs>
{
    enum Flag {NotYetInQueue=0, InQueue=1, NoLongerInQueue=2};
    enum Flag2 {Original, Added, New};
    enum ProcessedFlag {NotProcessed, Processed};

	Flag flag = NotYetInQueue;
	Flag2 flag2 = Original;
	ProcessedFlag processedFlag = NotProcessed;
public:
    MyHalfedge(){}

	inline void resetState()
	{
		flag = NotYetInQueue;
		flag2 = Original;
		processedFlag = NotProcessed;
	}

        /* Flag 1 */

	inline void setInQueue()
	{
	  flag=InQueue;
	}

	inline void removeFromQueue()
	{
	  assert(flag==InQueue);
	  flag=NoLongerInQueue;
	}

	/* Processed flag */

	inline void resetProcessedFlag()
	{
	  processedFlag = NotProcessed;
	}

	inline void setProcessed()
	{
		processedFlag = Processed;
	}

	inline bool isProcessed() const
	{
		return (processedFlag == Processed);
	}

	/* Flag 2 */

	inline void setAdded()
	{
	  assert(flag2 == Original);
	  flag2=Added;
	}

	inline void setNew()
	{
		assert(flag2 == Original);
		flag2 = New;
	}

	inline bool isAdded() const
	{
	  return flag2==Added;
	}

	inline bool isOriginal() const
	{
	  return flag2==Original;
	}

	inline bool isNew() const
	{
	  return flag2 == New;
	}
};

inline float triangle_area(const Point &p1, const Point &p2, const Point &p3){
	float x1 = p2.x()-p1.x();
	float y1 = p2.y()-p1.y();
	float z1 = p2.z()-p1.z();

	float x2 = p3.x()-p1.x();
	float y2 = p3.y()-p1.y();
	float z2 = p3.z()-p1.z();

	float x3 = y1*z2 - y2*z1;
	float y3 = z1*x2 - z2*x1;
	float z3 = x1*y2 - x2*y1;

	return sqrt(x3*x3 + y3*y3 + z3*z3)/2;
}

inline float triangle_area(const Triangle &tri){
	return triangle_area(tri[0], tri[1], tri[2]);
}

class replacing_group;

class MyTriangle;

// My face type has a vertex flag
template <class Refs>
class MyFace : public CGAL::HalfedgeDS_face_base<Refs>
{
    enum Flag {Unknown=0, Splittable=1, Unsplittable=2};
    enum ProcessedFlag {NotProcessed, Processed};

	Flag flag = Unknown;
	ProcessedFlag processedFlag = NotProcessed;

	Point removedVertexPos;
	float proxy_hausdorff_distance = 0.0;
	float hausdorff_distance = 0.0;
public:
    MyFace(){}

	inline void resetState()
	{
          flag = Unknown;
          processedFlag = NotProcessed;
	}

	inline void resetProcessedFlag()
	{
	  processedFlag = NotProcessed;
	}

	inline bool isConquered() const
	{
	  return (flag==Splittable ||flag==Unsplittable) ;
	}

	inline bool isSplittable() const
	{
	  return (flag==Splittable) ;
	}

	inline bool isUnsplittable() const
	{
	  return (flag==Unsplittable) ;
	}

	inline void setSplittable()
	{
	  assert(flag == Unknown);
	  flag=Splittable;
	}

	inline void setUnsplittable()
	{
	  assert(flag == Unknown);
	  flag=Unsplittable;
	}

	inline void setProcessedFlag()
	{
		processedFlag = Processed;
	}

	inline bool isProcessed() const
	{
		return (processedFlag == Processed);
	}

	inline Point getRemovedVertexPos() const
	{
		return removedVertexPos;
	}

	inline void setRemovedVertexPos(Point p)
	{
		removedVertexPos = p;
	}

public:

	inline pair<float, float> getHausdorfDistance(){
		return pair<float, float>(proxy_hausdorff_distance, hausdorff_distance);
	}

	inline void resetHausdorff(){
		hausdorff_distance = 0.0;
		proxy_hausdorff_distance = 0.0;
	}

	inline void setProxyHausdorff(float prh){
		proxy_hausdorff_distance = prh;
	}

	inline void setHausdorff(float hd){
		hausdorff_distance = hd;
	}

	inline void updateProxyHausdorff(float prh){
		proxy_hausdorff_distance = max(prh, proxy_hausdorff_distance);
	}

	inline void updateHausdorff(float hd){
		hausdorff_distance = max(hd, hausdorff_distance);
	}

	inline float getProxyHausdorff(){
		return proxy_hausdorff_distance;
	}

	inline float getHausdorff(){
		return hausdorff_distance;
	}
	replacing_group *rg = NULL;
	vector<Triangle> triangles;
	MyTriangle *tri = NULL;
};


struct MyItems : public CGAL::Polyhedron_items_3
{
    template <class Refs, class Traits>
    struct Face_wrapper {
        typedef MyFace<Refs> Face;
    };

	template <class Refs, class Traits>
    struct Vertex_wrapper {
        typedef MyVertex<Refs> Vertex;
    };

	template <class Refs, class Traits>
    struct Halfedge_wrapper {
        typedef MyHalfedge<Refs> Halfedge;
    };
};

// four types of hausdorff computing methods
enum Hausdorff_Computing_Type{
	HCT_NULL = 0,
	HCT_BVHTREE = 1,
	HCT_ASSOCIATE = 2,
	HCT_ASSOCIATE_CYLINDER = 3
};

class HiMesh: public CGAL::Polyhedron_3< MyKernel, MyItems >
{
	// Gate queues
	std::queue<Halfedge_handle> gateQueue;

	// Processing mode: 0 for compression and 1 for decompression.
	int i_mode;
	bool b_jobCompleted = false; // True if the job has been completed.

	unsigned i_curDecimationId = 0;
	unsigned i_nbDecimations;
	unsigned i_decompPercentage = 0;
	// Number of vertices removed during current conquest.
	unsigned i_nbRemovedVertices;

	// The vertices of the edge that is the departure of the coding and decoding conquests.
	Vertex_handle vh_departureConquest[2];
	// Geometry symbol list.
	std::deque<std::deque<Point> > geometrySym;

	std::deque<std::deque<unsigned char>> hausdorfSym;
	std::deque<std::deque<unsigned char>> proxyhausdorfSym;

	std::deque<std::deque<float>> hausdorfSym_float;
	std::deque<std::deque<float>> proxyhausdorfSym_float;

	// Connectivity symbol list.
	std::deque<std::deque<unsigned> > connectFaceSym;
	std::deque<std::deque<unsigned> > connectEdgeSym;


	// The compressed data;
	char *p_data;
	size_t dataOffset = 0; // the offset to read and write.

	aab mbb; // the bounding box
	SegTree *segment_tree = NULL;
	TriangleTree *triangle_tree = NULL;
	list<Segment> segments;
	list<Triangle> aabb_triangles;

	vector<MyTriangle *> original_facets;

	// Store the maximum Hausdorf Distance
	vector<pair<float, float>> globalHausdorfDistance;
	vector<Point> removedPoints;

	//
	unordered_set<replacing_group *> map_group;

	bool own_data = true;
public:

	int id = 0;
	// constructor for encoding
	HiMesh(string &str, bool completeop = false);

	// constructors for decoding
	HiMesh(char *data, size_t dsize, bool own_data = true);
	HiMesh(HiMesh *mesh): HiMesh(mesh->p_data, mesh->dataOffset, true){}
	~HiMesh();

	void encode(int lod = 0);
	void decode(int lod = 100);

	// Compression
	void startNextCompresssionOp();
	void RemovedVertexCodingStep();
	void InsertedEdgeCodingStep();
	void HausdorffCodingStep();

	void merge(unordered_set<replacing_group *> &reps, replacing_group *);
	Halfedge_handle vertexCut(Halfedge_handle startH);
	void encodeInsertedEdges(unsigned i_operationId);
	void encodeRemovedVertices(unsigned i_operationId);
	void encodeHausdorff(unsigned i_operationId);


	// Compression geometry and connectivity tests.
	bool isRemovable(Vertex_handle v) const;
	bool isConvex(const std::vector<Vertex_const_handle> & polygon) const;
	bool isPlanar(const std::vector<Vertex_const_handle> &polygon, float epsilon) const;
	bool willViolateManifold(const std::vector<Halfedge_const_handle> &polygon) const;
	float removalError(Vertex_const_handle v, const std::vector<Vertex_const_handle> &polygon) const;

	// Decompression
	void startNextDecompresssionOp();
	void RemovedVerticesDecodingStep();
	void InsertedEdgeDecodingStep();
	void HausdorffDecodingStep();
	void insertRemovedVertices();
	void removeInsertedEdges();

	// Utils
	Vector computeNormal(Halfedge_const_handle heh_gate) const;
	Vector computeNormal(const std::vector<Vertex_const_handle> & polygon) const;
	Vector computeNormal(Point p[3]) const;
	Vector computeNormal(Facet_const_handle f) const;
	Vector computeVertexNormal(Halfedge_const_handle heh) const;
	Point barycenter(Facet_const_handle f) const;
	Point barycenter(Halfedge_handle heh_gate) const;
	Point barycenter(const std::vector<Vertex_const_handle> &polygon) const;
	unsigned vertexDegreeNotNew(Vertex_const_handle vh) const;
	float triangleSurface(const Point p[]) const;
	float edgeLen(Halfedge_const_handle heh) const;
	float facePerimeter(const Face_handle fh) const;
	float faceSurface(Halfedge_handle heh) const;
	void pushHehInit();

	void alphaFolding();

	// IOs
	void writeFloat(float f);
	float readFloat();
	void writeInt16(int16_t i);
	int16_t readInt16();
	void writeuInt16(uint16_t i);
	uint16_t readuInt16();
	void writeInt(int i);
	int readInt();
	unsigned char readChar();
	void writeChar(unsigned char ch);
	void writePoint(Point &p);
	Point readPoint();

	void writeBaseMesh();
	void readBaseMesh();

	//3dpro
	void computeHausdorfDistance();
	void updateVFMap();
	bool isProtruding(const std::vector<Halfedge_const_handle> &polygon) const;
	void profileProtruding();

	inline aab get_mbb(){
		return mbb;
	}
	void compute_mbb();

	Polyhedron *to_polyhedron();
	Polyhedron *to_triangulated_polyhedron();
	Skeleton *extract_skeleton();
	vector<Point> get_skeleton_points(int num_skeleton_points);
	vector<Voxel *> generate_voxels_skeleton(int voxel_size);
	vector<Voxel *> voxelization(int voxel_size);
	string to_wkt();
	string to_off();
	void write_to_off(const char *path);
	void write_to_wkt(const char *path);

	float get_volume();
	size_t size_of_triangles();

	// get the elements
	size_t fill_segments(float *&segments);
	size_t fill_triangles(float *&triangles);
	size_t fill_hausdorf_distances(float *&hausdorf);
	pair<float, float> get_triangle_hausdorf(int tri_id);
	size_t fill_voxels(vector<Voxel *> &voxels);
	size_t fill_vertices(float *&vertices);
	list<Segment> get_segments();
	list<Triangle> get_triangles();
	map<float, Face_iterator> get_fits();
	list<Point> get_vertices();

	// query
	SegTree *get_aabb_tree_segment();
	TriangleTree *get_aabb_tree_triangle();
	Tree *get_aabb_tree();
	float distance(HiMesh *target);
	float distance_tree(HiMesh *target);
	float distance_tree(const Point &p);
	bool intersect(HiMesh *target);
	bool intersect_tree(HiMesh *target);

	void clear_aabb_tree();

	size_t size_of_edges();

	float getProxyHausdorffDistance();
	float getHausdorffDistance();

//	pair<float, float> getHausdorfDistance();
//	pair<float, float> getNextHausdorfDistance();

	bool is_compression_mode(){
		return i_mode == COMPRESSION_MODE_ID;
	}
	size_t get_data_size(){
		return dataOffset;
	}
	const char *get_data(){
		return p_data;
	}

	aab shift(float x_sft, float y_sft, float z_sft);
	aab shrink(float ratio);
	HiMesh *clone_mesh();

	map<Point, vector<MyTriangle *>> VFmap;

	// equals the number of points sampled for each triangle
	static uint sampling_rate;

	static int calculate_method;
	static bool use_byte_coding;

	void sample_points(const HiMesh::Face_iterator &fit, unordered_set<Point> &points, float area_unit);
	void sample_points(const Triangle &tri, unordered_set<Point> &points, float area_unit);
	void sample_points(unordered_set<Point> &points);
	inline float get_sample_density(){
		uint num_triangle = this->size_of_triangles();
		return area()/(num_triangle*sampling_rate);
	}
	float area();
};

class MyTriangle{
public:
	MyTriangle(const Triangle &t, int i=0){
		tri = t;
		processed = false;
		id = i;
	}
	~MyTriangle(){
		sampled_points.clear();
		facets.clear();
	}
	// owned
	int id;
	Triangle tri;
	unordered_set<Point> sampled_points;

	bool processed = false;
	// associated
	vector<HiMesh::Face_iterator> facets;
	void add_facet(HiMesh::Face_iterator fit){
		bool contained = false;
		for(HiMesh::Face_iterator f:facets){
			if(f==fit){
				return;
			}
		}
		facets.push_back(fit);
	}
	void reset(){
		processed = false;
		facets.clear();
	}
};

class replacing_group{
public:
	replacing_group(){
		//cout<<this<<" is constructed"<<endl;
		id = counter++;
		alive++;
	}
	~replacing_group(){
		removed_vertices.clear();
		alive--;
	}

	void print(){
		log("%5d (%2d refs %4d alive) - removed_vertices: %ld", id, ref, alive, removed_vertices.size());
	}

	unordered_set<Point> removed_vertices;
	//unordered_set<Triangle> removed_triangles;
	int id;
	int ref = 0;

	static int counter;
	static int alive;
};

enum Decoding_Type{
	COMPRESSED,
	RAW
};

/*
 * a wrapper for describing a mesh. for some cases we may not need to truly
 * parse the mesh out from disk, but use the bounding box is good enough.
 * */
class HiMesh_Wrapper{

	map<int, float> hausdorffs;
	map<int, float> proxyhausdorffs;
public:
	Decoding_Type type;
	// the data mode

	// description of the polyhedron
	size_t id = -1;
	weighted_aab box;
	vector<Voxel *> voxels;

	// used for retrieving compressed data from disk
	HiMesh *mesh = NULL;
	char *data_buffer = NULL;
	size_t data_size = 0;

	char *meta_buffer = NULL;
	size_t meta_size = 0;

	pthread_mutex_t lock;
	vector<HiMesh_Wrapper *> results;
	int cur_lod = -1;

public:
	HiMesh_Wrapper(char *dt, size_t id, Decoding_Type t = COMPRESSED);
	~HiMesh_Wrapper();

	void decode_to(int lod);
	void clear_voxels();

	// for RAW data mode
	size_t get_voxel_offset(int id, int lod);
	size_t get_voxel_size(int id, int lod);

	float getHausdorffDistance();
	float getProxyHausdorffDistance();

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
void write_triangles(vector<Triangle> &triangles, const char *path);
void write_triangles(vector<Triangle *> &triangles, const char *path);
string read_off_stdin();
string polyhedron_to_wkt(Polyhedron *poly);

// some utility functions to operate mesh polyhedrons
extern HiMesh *read_mesh(bool complete_compress = false);
extern HiMesh *read_mesh(char *path, bool complete_compress = false);
extern vector<HiMesh *> read_meshes(const char *path, size_t maxnum = LONG_MAX);

extern HiMesh *poly_to_mesh(Polyhedron *poly);

float get_volume(Polyhedron *polyhedron);

Polyhedron *parse_polyhedron(string &str);
Polyhedron *read_polyhedron();
extern Polyhedron *read_polyhedron(const char *path);
vector<Polyhedron *> read_polyhedrons(const char *path, size_t load_num = LONG_MAX);

// get the Euclid distance of two points
inline float get_distance(const Point &p1, const Point &p2){
	float dist = 0;
	for(int i=0;i<3;i++){
		dist += (p2[i]-p1[i])*(p2[i]-p1[i]);
	}
	return dist;
}

}


#endif
