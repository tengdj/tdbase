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
	data_buffer = dt;
	data_size = sz;
	size_t vnum = *(size_t *)(data_buffer + data_size);
	voxel_meta_buffer = data_buffer+data_size;
	voxel_meta_size = sizeof(vnum);
	if(type == COMPRESSED){
		// the mesh reuses the memory space stored in the Tile class
		mesh = new HiMesh(data_buffer, data_size, false);
		mesh->id = id;
		for(int i=0;i<vnum;i++){
			Voxel *v = new Voxel();
			memcpy(v->low, voxel_meta_buffer+voxel_meta_size, 3*sizeof(float));
			voxel_meta_size += 3*sizeof(float);
			memcpy(v->high, voxel_meta_buffer+voxel_meta_size, 3*sizeof(float));
			voxel_meta_size += 3*sizeof(float);
			memcpy(v->core, voxel_meta_buffer+voxel_meta_size, 3*sizeof(float));
			voxel_meta_size += 3*sizeof(float);
			voxels.push_back(v);
			box.update(*v);
		}
	}else{
		// record the decoded data for different voxels in different LODs
		for(int i=0;i<vnum;i++){
			Voxel *v = new Voxel();
			// load the space information
			memcpy(v->low, voxel_meta_buffer+voxel_meta_size, 3*sizeof(float));
			voxel_meta_size += 3*sizeof(float);
			memcpy(v->high, voxel_meta_buffer+voxel_meta_size, 3*sizeof(float));
			voxel_meta_size += 3*sizeof(float);
			memcpy(v->core, voxel_meta_buffer+voxel_meta_size, 3*sizeof(float));
			voxel_meta_size += 3*sizeof(float);

			// load the offset and volume information for varying LODs
			for(int lod=0;lod<=100;lod+=10){
				size_t of = *(size_t *)(voxel_meta_buffer+voxel_meta_size);
				voxel_meta_size += sizeof(size_t);
				size_t vl = *(size_t *)(voxel_meta_buffer+voxel_meta_size);
				voxel_meta_size += sizeof(size_t);
				v->offset_lod[lod] = of;
				v->volumn_lod[lod] = vl;

				//printf("%ld\t%ld\t%ld\t%ld\t%ld\n", id, i, lod, v->offset_lod[lod], v->volumn_lod[lod]);
			}

			voxels.push_back(v);
			box.update(*v);
		}
	}

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

void HiMesh_Wrapper::decode_to(int lod){
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
	}else{
		// for RAW data mode, simply link pointers instead of do the decoding job
		// as data already been stored a
		for(int i=0;i<voxels.size();i++){
			size_t offset = get_voxel_offset(i, lod);
			size_t size = get_voxel_size(i, lod);
			voxels[i]->external_load((float *)(data_buffer+offset), (float *)(data_buffer+offset+9*size*sizeof(float)), size);
		}
	}
	if(global_ctx.verbose>=2){
		for(int i=0;i<voxels.size();i++){
			printf("%d %d %d %d %d\n", id, i, lod, voxels[i]->offset_lod[lod], voxels[i]->volumn_lod[lod]);
			voxels[i]->print();
		}
	}
}

size_t HiMesh_Wrapper::get_voxel_offset(int id, int lod){
	return voxels[id]->offset_lod[lod];
}
size_t HiMesh_Wrapper::get_voxel_size(int id, int lod){
	return voxels[id]->volumn_lod[lod];
}

void HiMesh_Wrapper::report_result(HiMesh_Wrapper *result){
	pthread_mutex_lock(&lock);
	results.push_back(result);
	pthread_mutex_unlock(&lock);
	//log("%d find %d", id, result->id);
}

}
