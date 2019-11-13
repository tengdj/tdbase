/*
 * spatial.h
 *
 * any CGAL spatial related implementations
 * will include this header
 *
 *  Created on: Nov 12, 2019
 *      Author: teng
 */

#ifndef SPATIAL_CGAL_H_
#define SPATIAL_CGAL_H_

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/intersections.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/convex_decomposition_3.h>
#include <CGAL/Tetrahedron_3.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/extract_mean_curvature_flow_skeleton.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/boost/graph/split_graph_into_polylines.h>
#include <CGAL/OFF_to_nef_3.h>

#include <boost/foreach.hpp>

#include "../PPMC/ppmc.h"

using namespace CGAL;

// some local definition

// MyKernel is inherited from the PPMC project
typedef CGAL::Polyhedron_3<MyKernel>	Polyhedron;
typedef Polyhedron::Halfedge_handle Halfedge_handle;
typedef CGAL::Surface_mesh<Point>                             Triangle_mesh;
typedef boost::graph_traits<Triangle_mesh>::vertex_descriptor vertex_descriptor;
typedef CGAL::Mean_curvature_flow_skeletonization<Triangle_mesh> Skeletonization;
typedef Skeletonization::Skeleton                             Skeleton;
typedef Skeleton::vertex_descriptor                           Skeleton_vertex;
typedef Skeleton::edge_descriptor                             Skeleton_edge;
typedef MyMesh::Vertex_iterator Vertex_iterator;

namespace hispeed{

class mbb{
public:
	float min[3];
	float max[3];
	mbb(){
		for(int i=0;i<3;i++){
			min[i] = DBL_MAX;
			max[i] = -DBL_MAX;
		}
	}
	mbb(Point min, Point max){
		for(int i=0;i<3;i++){
			this->min[i] = min[i];
			this->max[i] = max[i];
		}
	}
	mbb(std::vector<Point> &points){
		for(Point p:points){
			update(p);
		}
	}

	void update(Point &p){
		for(int i=0;i<3;i++){
			if(min[i]>p[i]){
				min[i] = p[i];
			}
			if(max[i]<p[i]){
				max[i] = p[i];
			}
		}
	}

	friend std::ostream&
	operator<<(std::ostream& os, const mbb &p){
		for(int i=0;i<3;i++){
			os<<p.min[i]<<" ";
		}
		os<<"-> ";
		for(int i=0;i<3;i++){
			os<<p.max[i]<<" ";
		}
		return os;
	}


};



Skeleton *extract_skeleton(MyMesh *current_mesh);
void get_skeleton_points(Skeleton &skeleton, std::vector<Point> &P);
void get_skeleton_edges(Skeleton &skeleton);

Polyhedron *make_cube(mbb box);

// some utility functions to operate mesh polyhedrons
void print_mesh(Polyhedron *mesh);
void print_mesh_file(Polyhedron *mesh, char *path);
MyMesh *read_mesh();
MyMesh *decompress_mesh(MyMesh *compressed, int lod);

}



#endif /* SPATIAL_CGAL_H_ */
