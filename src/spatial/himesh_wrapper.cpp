/*
 * himesh_wrapper.cpp
 *
 *  Created on: Jun 24, 2022
 *      Author: teng
 */


#include "himesh.h"

namespace tdbase{

/*
 * himesh wrapper functions
 * */

HiMesh_Wrapper::HiMesh_Wrapper(char *dt, size_t i, Decoding_Type t){
	type = t;
	id = i;
	box.id = i;
	data_size = *(size_t *)(dt);

	data_buffer = dt + sizeof(size_t);
	meta_buffer = data_buffer + data_size;

	size_t vnum = *(size_t *)meta_buffer;
	meta_size = sizeof(vnum);
	if(type == COMPRESSED){
		// the mesh reuses the memory space stored in the Tile class
		mesh = new HiMesh(data_buffer, data_size, false);

		for(int i=0;i<vnum;i++){
			Voxel *v = new Voxel();
			memcpy(v->low, meta_buffer+meta_size, 3*sizeof(float));
			meta_size += 3*sizeof(float);
			memcpy(v->high, meta_buffer+meta_size, 3*sizeof(float));
			meta_size += 3*sizeof(float);
			memcpy(v->core, meta_buffer+meta_size, 3*sizeof(float));
			meta_size += 3*sizeof(float);
			voxels.push_back(v);
			box.update(*v);
		}
	}else{
		for(int lod=20;lod<=100;lod+=20){
			this->hausdorffs[lod] = *(float *)(meta_buffer+meta_size);
			meta_size += sizeof(float);
			this->proxyhausdorffs[lod] = *(float *)(meta_buffer+meta_size);
			meta_size += sizeof(float);
		}

		// record the decoded data for different voxels in different LODs
		for(int i=0;i<vnum;i++){
			Voxel *v = new Voxel();
			// load the space information
			memcpy(v->low, meta_buffer+meta_size, 3*sizeof(float));
			meta_size += 3*sizeof(float);
			memcpy(v->high, meta_buffer+meta_size, 3*sizeof(float));
			meta_size += 3*sizeof(float);
			memcpy(v->core, meta_buffer+meta_size, 3*sizeof(float));
			meta_size += 3*sizeof(float);

			// load the offset and volume information for varying LODs
			for(int lod=20;lod<=100;lod+=20){
				size_t of = *(size_t *)(meta_buffer+meta_size);
				meta_size += sizeof(size_t);
				size_t vl = *(size_t *)(meta_buffer+meta_size);
				meta_size += sizeof(size_t);
				v->offset_lod[lod] = of;
				v->volume_lod[lod] = vl;
				//printf("%ld\t%ld\t%ld\t%ld\t%ld\n", id, i, lod, v->offset_lod[lod], v->volumn_lod[lod]);
			}

			voxels.push_back(v);
			box.update(*v);
		}
	}
	pthread_mutex_init(&lock, NULL);
}

// manually
HiMesh_Wrapper::HiMesh_Wrapper(map<int, HiMesh *> &ms){
	type = MULTIMESH;
	for(auto a:ms){
		meshes[a.first] = a.second;
	}
	assert(ms.find(100)!=ms.end());

	voxels = ms[100]->generate_voxels_skeleton();
	for(Voxel *v:voxels){
		box.update(*v);
	}

}

HiMesh_Wrapper::HiMesh_Wrapper(HiMesh *m){
	type = COMPRESSED;
	mesh = m;
	voxels = m->generate_voxels_skeleton();
	m->encode();
}

HiMesh_Wrapper::~HiMesh_Wrapper(){
	for(Voxel *v:voxels){
		delete v;
	}
	voxels.clear();
	if(mesh){
		delete mesh;
	}
	for(auto a:meshes){
		delete a.second;
	}
	meshes.clear();
	results.clear();
}

void HiMesh_Wrapper::decode_to(int lod){
	if(lod <= cur_lod){
		return;
	}
	cur_lod = lod;
	for(Voxel *v:voxels){
		v->clear();
	}

	if(type == MULTIMESH){
		assert(meshes.find(cur_lod)!=meshes.end());
		get_mesh()->fill_voxels(voxels);
	}else if(type == COMPRESSED){
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
	if(global_ctx.verbose>=3){
		for(int i=0;i<voxels.size();i++){
			printf("decode_to: id: %ld\t voxel_id: %d\t lod: %d\t offset: %ld\t volume: %ld\n", id, i, lod, voxels[i]->offset_lod[lod], voxels[i]->volume_lod[lod]);
			voxels[i]->print();
		}
	}
}

float HiMesh_Wrapper::getHausdorffDistance(){
	if(type == COMPRESSED){
		return mesh->getHausdorffDistance();
	}else{
		assert(hausdorffs.find(cur_lod)!=hausdorffs.end());
		return hausdorffs[cur_lod];
	}
}
float HiMesh_Wrapper::getProxyHausdorffDistance(){
	if(type == COMPRESSED){
		return mesh->getProxyHausdorffDistance();
	}else{
		assert(hausdorffs.find(cur_lod)!=hausdorffs.end());
		return proxyhausdorffs[cur_lod];
	}
}

size_t HiMesh_Wrapper::get_voxel_offset(int id, int lod){
	return voxels[id]->offset_lod[lod];
}
size_t HiMesh_Wrapper::get_voxel_size(int id, int lod){
	return voxels[id]->volume_lod[lod];
}

void HiMesh_Wrapper::report_result(HiMesh_Wrapper *result){
	pthread_mutex_lock(&lock);
	results.push_back(result);
	if (global_ctx.print_result) {
		printf("%ld %ld\n", id, result->id);
	}
	pthread_mutex_unlock(&lock);
}

}
