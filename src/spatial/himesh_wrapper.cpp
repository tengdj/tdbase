/*
 * himesh_wrapper.cpp
 *
 *  Created on: Jun 24, 2022
 *      Author: teng
 */


#include "himesh.h"

namespace hispeed{

/*
 * himesh wrapper functions
 * */

HiMesh_Wrapper::HiMesh_Wrapper(char *dt, size_t sz, size_t i, Decoding_Type t){
	type = t;
	id = i;
	box.id = i;
	// the mesh use the memory space stored in the Tile class, avoid one copy
	mesh = new HiMesh(dt, sz, false);
	mesh->id = id;
	data_buffer = dt;
	data_size = sz;
	pthread_mutex_init(&lock, NULL);
}

HiMesh_Wrapper::~HiMesh_Wrapper(){
	for(Voxel *v:voxels){
		delete v;
	}
	voxels.clear();
	if(mesh){
		delete mesh;
	}
	results.clear();
}

void HiMesh_Wrapper::decode_to(uint lod){
	if(lod <= cur_lod){
		return;
	}
	for(Voxel *v:voxels){
		v->clear();
	}
	cur_lod = lod;
	if(type == COMPRESSED){
		assert(mesh);
		mesh->decode(lod);
		mesh->fill_voxels(voxels);
	}
}

void HiMesh_Wrapper::disable_innerpart(){
	if(voxels.size()>1){
		for(Voxel *v:voxels){
			delete v;
		}
		voxels.clear();
		Voxel *v = new Voxel();
		v->set_box(box);
		for(int i=0;i<3;i++){
			v->core[i] = (v->high[i]-v->low[i])/2+v->low[i];
		}
		voxels.push_back(v);
	}
}

void HiMesh_Wrapper::report_result(HiMesh_Wrapper *result){
	pthread_mutex_lock(&lock);
	results.push_back(result);
	pthread_mutex_unlock(&lock);
	//log("%d find %d", id, result->id);
}

}
