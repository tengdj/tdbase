/*
 * skeleton.h
 *
 *  Created on: Nov 12, 2019
 *      Author: teng
 */

#ifndef SKELETON_H_
#define SKELETON_H_

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
typedef CGAL::Surface_mesh<Point>                             Triangle_mesh;
typedef boost::graph_traits<Triangle_mesh>::vertex_descriptor vertex_descriptor;
typedef CGAL::Mean_curvature_flow_skeletonization<Triangle_mesh> Skeletonization;
typedef Skeletonization::Skeleton                             Skeleton;
typedef Skeleton::vertex_descriptor                           Skeleton_vertex;
typedef Skeleton::edge_descriptor                             Skeleton_edge;
namespace hispeed{


Skeleton *extract_skeleton(MyMesh *current_mesh);
void get_skeleton_points(Skeleton &skeleton, std::vector<Point> &P);
void get_skeleton_edges(Skeleton &skeleton);

}



#endif /* SKELETON_H_ */
