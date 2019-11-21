/*
 * voxel_group.cpp
 *
 *  Created on: Nov 20, 2019
 *      Author: teng
 */



#include "tile.h"


namespace hispeed{


Voxel_group::Voxel_group(){
	for(int i=0;i<10;i++){
		decoded_data[i] = NULL;
	}
};

Voxel_group::~Voxel_group(){
	// release the on-demand decoded data
	release_data();
	for(Voxel *v:voxels){
		assert(v);
		delete v;
	}
	voxels.clear();
}

void Voxel_group::add_voxel(aab box){
	Voxel *v = new Voxel();
	v->box = box;
	v->id = voxels.size();
	voxels.push_back(v);
	v->group = this;
	this->box.update(box);
}

bool Voxel_group::persist(FILE *fs){
	fwrite(&offset, 1, sizeof(long), fs);
	fwrite(&data_size, 1, sizeof(long), fs);
	uint size = voxels.size();
	fwrite(&size, 1, sizeof(uint), fs);
	for(Voxel *v:voxels){
		fwrite(v->box.min, 3, sizeof(float), fs);
		fwrite(v->box.max, 3, sizeof(float), fs);
	}
	return true;
}

bool Voxel_group::load(FILE *fs){
	fread((void *)&offset, sizeof(long), 1, fs);
	fread((void *)&data_size, sizeof(long), 1, fs);
	uint size = 0;
	fread((void *)&size, sizeof(uint), 1, fs);
	// load the voxels
	for(uint i=0;i<size;i++){
		aab tmpb;
		fread((void *)tmpb.min, sizeof(float), 3, fs);
		fread((void *)tmpb.max, sizeof(float), 3, fs);
		add_voxel(tmpb);
	}
	return true;
}



void Voxel_group::decode(int lod){

	// load the mesh object from disk
	if(mesh == NULL){
		tile->retrieve_mesh(id);
	}
	assert(mesh);
	assert(decoded_data[lod] = NULL);
	decoded_data[lod] = mesh->decode_lod(lod);

	// set the data pointers to the
	// related voxels
	float *cur_data = decoded_data[lod];
	for(int i=0;i<voxels.size();i++){
		voxels[i]->set_data(lod, cur_data);
		// the first element contains the size
		cur_data += 1+(int)(*cur_data);
	}
}

void Voxel_group::release_data(){
	for(Voxel *v:voxels){
		v->clear();
	}
	if(mesh!=NULL){
		delete mesh;
		mesh = NULL;
		delete mesh_data;
		mesh_data = NULL;
		for(int i=0;i<10;i++){
			if(decoded_data[i]){
				delete decoded_data[i];
				decoded_data[i] = NULL;
			}
		}
	}
}


void Voxel_group::print(){
	cout<<"\tvoxel group "<<id<<endl;
	for(Voxel *v:voxels){
		v->print();
	}
}

}
