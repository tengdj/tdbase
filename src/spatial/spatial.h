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

#include "../PPMC/ppmc.h"
#include "../util/util.h"

using namespace CGAL;
using namespace std;

// some local definition

namespace hispeed{

Polyhedron *make_cube(mbb box);

void write_polyhedron(Polyhedron *mesh, char *path);
// some utility functions to operate mesh polyhedrons
MyMesh *read_mesh();
MyMesh *decompress_mesh(MyMesh *compressed, int lod);

}



#endif /* SPATIAL_CGAL_H_ */
