/*
 * tile.cpp
 *
 *  Created on: Nov 14, 2019
 *      Author: teng
 *
 *  the implementation of functions manipulating
 *  disk files and memory space which stores the
 *  compressed polyhedra and light weight index
 *  of each polyhedra. Each file contains the
 *  information for one single tile
 *
 */


#include "tile.h"


namespace hispeed{


void Voxel_group::decode(int lod){

	// load the mesh object from disk
	if(mesh == NULL){

	}

	float *decoded_data;
	mesh->decode_lod(lod, &decoded_data);

	// set the data pointers to the
	// related voxels
	for(int i=0;i<voxels.size();i++){
		int size = (int)*decoded_data;
		voxels[i]->set_data(lod, decoded_data);
		decoded_data += 1+size;
	}

}

void Voxel_group::release_data(){
	if(mesh!=NULL){
		delete mesh;
		mesh = NULL;
	}
	for(Voxel *v:voxels){
		v->clear();
	}
	this->data_size = 0;
}


/*
 * persist the information of the tile
 * to disk
 *
 * */
bool Tile::persist(){
	if(fs!=NULL){
		fclose(fs);
	}
	fs = fopen(file_path.c_str(), "wb+");
	if(fs==NULL){
		cerr<<"failed to open file "<<file_path<<endl;
		return false;
	}

	fwrite(&id, 1, sizeof(long), fs);
	fwrite(&space.min, 3, sizeof(float), fs);
	fwrite(&space.max, 3, sizeof(float), fs);
	uint size = this->voxels;
	fwrite(&size, 1, size*sizeof(uint), fs);
	for(Voxel v:voxels){

	}
	fclose(fs);
	fs = NULL;
	return true;
}

/*
 * load the persisted information from disk
 * to construct the tile. Note that the information
 * of the voxels in the tile are parsed lazily only
 * on-demand. Thus the load function only load the description
 * info.
 * */
bool Tile::load(){
	if(fs==NULL){
		fs = fopen(file_path.c_str(), "r");
		if(fs==NULL){
			cerr<<"failed to open file "<<file_path<<endl;
			return false;
		}
	}


}



}

