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

// load meta data from file
// and construct the hierarchy structure
// tile->mesh->voxels->triangle/edges
Tile::Tile(std::string path){
	if(!hispeed::file_exist(path.c_str())){
		cerr<<path<<" does not exist"<<endl;
		exit(-1);
	}
	dt_fs = fopen(path.c_str(), "r");
	if(!dt_fs){
		cerr<<path<<" can not be opened"<<endl;
		exit(-1);
	}
	load();
	cout<<objects.size()<<" polyhedra in this tile"<<endl;
}

Tile::~Tile(){
	for(HiMesh_Wrapper *h:objects){
		delete h;
	}
	// close the data file pointer
	if(dt_fs!=NULL){
		fclose(dt_fs);
		dt_fs = NULL;
	}
}

// parse the metadata
bool Tile::load(){
	assert(dt_fs);
	size_t dsize = 0;
	long offset = 0;
	int index = 0;
	fseek(dt_fs, 0, SEEK_SET);
	while(fread((void *)&dsize, sizeof(size_t), 1, dt_fs)>0){
		fseek(dt_fs, dsize, SEEK_CUR);
		HiMesh_Wrapper *w = new HiMesh_Wrapper();
		offset += sizeof(size_t);
		w->offset = offset;
		w->data_size = dsize;
		w->id = index++;
		// read the voxels into the wrapper
		fread((void *)&dsize, sizeof(size_t), 1, dt_fs);
		for(int i=0;i<dsize;i++){
			Voxel *v = new Voxel();
			fread((void *)v->box.min, sizeof(float), 3, dt_fs);
			fread((void *)v->box.max, sizeof(float), 3, dt_fs);
			fread((void *)v->core, sizeof(float), 3, dt_fs);
			w->voxels.push_back(v);
			w->box.update(v->box);
		}
		objects.push_back(w);
		box.update(w->box);
		// update the offset for next
		offset += w->data_size+sizeof(size_t)+9*sizeof(float)*dsize;
	}
	return true;
}

// retrieve the mesh of the voxel group with ID id on demand
void Tile::retrieve_mesh(int id){
	assert(id>=0&&id<objects.size());
	HiMesh_Wrapper *wrapper = objects[id];
	assert(wrapper->mesh==NULL);
	char *mesh_data = new char[wrapper->data_size];
	assert(dt_fs);
	fseek(dt_fs, wrapper->offset, SEEK_SET);
	assert(wrapper->data_size==
			fread(mesh_data, sizeof(char), wrapper->data_size, dt_fs));
	wrapper->mesh = new HiMesh(mesh_data, wrapper->data_size);
	delete mesh_data;
}


void Tile::print(){
	int index = 0;
	for(HiMesh_Wrapper *w:objects){
		cout<<index++<<"\t"<<w->box<<endl;
	}
	cout<<endl;
}

}

