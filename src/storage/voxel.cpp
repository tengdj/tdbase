/*
 * voxel.cpp
 *
 *  Created on: Nov 20, 2019
 *      Author: teng
 */

#include "tile.h"


namespace hispeed{

void Voxel::clear(){
	for(int i=0;i<10;i++){
		data[i] = NULL;
		size[i] = 0;
	}
}
Voxel::Voxel(){
	clear();
};
Voxel::~Voxel(){
	clear();
}

int Voxel::get_size(int lod){
	return size[lod];
}
int Voxel::get_total_size(){
	int total_size = 0;
	for(int i=0;i<10;i++){
		total_size += size[i];
	}
	return total_size;
}

void Voxel::set_data(int lod, float *target_data){
	size[lod] = (int)(*target_data);
	data[lod] = target_data+1;
}

float *Voxel::get_data(int lod){
	assert(lod>=0&&lod<10);
	if(data[lod]==NULL){
		group->decode(lod);
	}
	return data[lod];
}

void Voxel::print(){
	cout<<"\t\t"<<box<<endl;
}

}


