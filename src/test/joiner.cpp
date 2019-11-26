/*
 * joiner.cpp
 *
 *  Created on: Nov 20, 2019
 *      Author: teng
 */



#include "../join/SpatialJoin.h"
#include "../storage/tile.h"
#include "../geometry/geometry.h"
using namespace std;
using namespace hispeed;

int main(int argc, char **argv){
	Tile *tile = new Tile("test");
	SpatialJoin *joiner = new SpatialJoin(tile, tile);
	joiner->formalize_computing();

	delete tile;
	delete joiner;

}
