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

	Tile *tile = new Tile("nuclei.dt");
	report_time("load tile", start);
	SpatialJoin *joiner = new SpatialJoin(tile, tile);
	joiner->formalize_computing();
	report_time("total join", start);

	delete tile;
	delete joiner;
}
