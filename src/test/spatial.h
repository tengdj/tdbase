/*
 * spatial.h
 *
 *  Created on: Oct 24, 2019
 *      Author: teng
 */

#ifndef SRC_TEST_SPATIAL_H_
#define SRC_TEST_SPATIAL_H_

#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <map>
#include <cstdlib>
#include <getopt.h>
#include <time.h>
#include <sys/shm.h>
#include <unordered_map>
#include <stdlib.h>

#include <boost/algorithm/string/replace.hpp>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/intersections.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/convex_decomposition_3.h>
#include <CGAL/Tetrahedron_3.h>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/extract_mean_curvature_flow_skeleton.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/boost/graph/split_graph_into_polylines.h>
#include <fstream>

#include <CGAL/Surface_mesh.h>


#include <boost/foreach.hpp>
#include <CGAL/Delaunay_triangulation_3.h>

#include <CGAL/OFF_to_nef_3.h>



//compressed data

#include <spatialindex/SpatialIndex.h>




#endif /* SRC_TEST_SPATIAL_H_ */
