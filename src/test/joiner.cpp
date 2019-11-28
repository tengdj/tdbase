/*
 * joiner.cpp
 *
 *  Created on: Nov 20, 2019
 *      Author: teng
 */



#include "../join/SpatialJoin.h"
#include "../storage/tile.h"
#include "../geometry/geometry.h"
#include <vector>
using namespace std;
using namespace hispeed;

int main(int argc, char **argv){
	struct timeval start = get_cur_time();
//	Tile *tile = new Tile(argv[1]);
//	SpatialJoin *joiner = new SpatialJoin(tile, tile);
//	joiner->formalize_computing();
//
//	delete tile;
//	delete joiner;
//	report_time("total join", start);
	Tile *tile = new Tile("vessel.dt");
	report_time("load tile", start);
	report_time("get mesh", start);
	HiMesh_Wrapper *wrapper = tile->get_mesh_wrapper(0);
	tile->get_mesh(0, 10);
	wrapper->fill_voxels(10);
	report_time("fill segments", start);

	delete tile;

}
