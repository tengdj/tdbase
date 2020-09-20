/*
 * towkt.cpp
 *
 *  Created on: Feb 15, 2020
 *      Author: teng
 */

#include "../spatial/spatial.h"
#include "../spatial/himesh.h"
#include "../storage/tile.h"

using namespace hispeed;

int main(int argc, char **argv){

	Tile *tile = new Tile(argv[1]);
	tile->retrieve_all();
	tile->advance_all(100);
	for(int i=0;i<tile->num_objects();i++){
		cout<<i<<"|"<<tile->get_mesh(i)->to_wkt()<<endl;
	}
	delete tile;
}


