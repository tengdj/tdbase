
#include "../spatial/himesh.h"
#include "../storage/tile.h"
#include "../PPMC/ppmc.h"
#include "../include/util.h"
#include <algorithm>
#include "zlib.h"

/*
 * adjust the size and the position of a polyhedron
 * */
int main(int argc, char **argv){
	Polyhedron *poly = hispeed::read_polyhedron();
	struct timeval start = get_cur_time();
	hispeed::cgal_simplification(poly, 0.1);
	logt("simplify", start);
	hispeed::write_polyhedron(poly, "/gisdata/simplified.off");
	delete poly;
}
