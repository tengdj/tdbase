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
#include <CGAL/OFF_to_nef_3.h>

using namespace CGAL;
using namespace std;

typedef CGAL::Nef_polyhedron_3<MyKernel> Nef_polyhedron;
typedef CGAL::Triangulation_3<MyKernel> Triangulation;
typedef Triangulation::Tetrahedron 	Tetrahedron;
typedef Nef_polyhedron::Volume_const_iterator Volume_const_iterator;

// some local definition

namespace hispeed{

Polyhedron *make_cube(aab box);

void write_polyhedron(Polyhedron *mesh, const char *path);
void write_polyhedron(Polyhedron *mesh, int id);
// some utility functions to operate mesh polyhedrons
extern MyMesh *get_mesh(string input, bool complete_compress = false);
extern MyMesh *read_mesh();
extern MyMesh *decompress_mesh(MyMesh *compressed, int lod);
extern MyMesh *decompress_mesh(char *data, size_t length, bool complete_operation = false);

inline void replace_bar(string &input){
	boost::replace_all(input, "|", "\n");
}

inline string read_polyhedron_str(){
	string input_line = read_line();
	boost::replace_all(input_line, "|", "\n");
	return input_line;
}
float get_volume(Polyhedron *polyhedron);
Polyhedron *read_polyhedron();

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
		dist += abs(p2[i]-p1[i]);
	}
	return dist;
}

}



#endif /* SPATIAL_CGAL_H_ */
