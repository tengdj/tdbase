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
#include "../geometry/aab.h"

using namespace CGAL;
using namespace std;

// some local definition

namespace hispeed{

Polyhedron *make_cube(aab box);

void write_polyhedron(Polyhedron *mesh, const char *path);
// some utility functions to operate mesh polyhedrons
extern MyMesh *get_mesh(string input, bool complete_compress = false);
MyMesh *read_mesh();
MyMesh *decompress_mesh(MyMesh *compressed, int lod);
extern MyMesh *compress_mesh(char *data, size_t length, bool complete_operation = false);


}



#endif /* SPATIAL_CGAL_H_ */
