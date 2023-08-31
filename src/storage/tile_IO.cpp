/*
 * tile_IO.cpp
 *
 *  Created on: Apr 4, 2023
 *      Author: teng
 */

#include "tile.h"

namespace hispeed{



//// parse the metadata
//bool Tile::parse_raw(){
//	FILE *dt_fs = fopen(tile_path.c_str(), "r");
//	assert(dt_fs);
//	size_t dsize = 0;
//	long offset = 0;
//	size_t index = 0;
//	fseek(dt_fs, 0, SEEK_SET);
//	while(fread((void *)&dsize, sizeof(size_t), 1, dt_fs)>0){
//		offset += sizeof(size_t);
//
//		// create a wrapper with the meta information
//		HiMesh_Wrapper *w = new HiMesh_Wrapper();
//		w->offset = offset;
//		w->data_size = dsize;
//		w->id = index++;
//		w->box.id = w->id;
//
//		// next dsize bytes are for the mesh, skip them
//		fseek(dt_fs, dsize, SEEK_CUR);
//		// read the voxels into the wrapper
//		fread((void *)&dsize, sizeof(size_t), 1, dt_fs);
//		for(int i=0;i<dsize;i++){
//			Voxel *v = new Voxel();
//			fread((void *)v->low, sizeof(float), 3, dt_fs);
//			fread((void *)v->high, sizeof(float), 3, dt_fs);
//			fread((void *)v->core, sizeof(float), 3, dt_fs);
//			w->voxels.push_back(v);
//			w->box.update(*v);
//		}
//		objects.push_back(w);
//		space.update(w->box);
//
//		// update the offset for next object
//		offset += w->data_size+sizeof(size_t)+9*sizeof(float)*dsize;
//	}
//	fclose(dt_fs);
//	return true;
//}
//
//// persist the meta data for current tile as a cache
//bool Tile::persist(string path){
//	FILE *mt_fs = fopen(path.c_str(), "wb+");
//	assert(mt_fs);
//	for(HiMesh_Wrapper *w:objects){
//		fwrite((void *)&w->offset, sizeof(size_t), 1, mt_fs);
//		fwrite((void *)&w->data_size, sizeof(size_t), 1, mt_fs);
//		size_t dsize = w->voxels.size();
//		fwrite((void *)&dsize, sizeof(size_t), 1, mt_fs);
//		for(Voxel *v:w->voxels){
//			fwrite((void *)v->low, sizeof(float), 3, mt_fs);
//			fwrite((void *)v->high, sizeof(float), 3, mt_fs);
//			fwrite((void *)v->core, sizeof(float), 3, mt_fs);
//		}
//	}
//	fclose(mt_fs);
//	return true;
//}
//
//// load from the cached meta data
//bool Tile::load(string path, int capacity){
//	FILE *mt_fs = fopen(path.c_str(), "r");
//	if(mt_fs==NULL){
//		log("%s cannot be opened, error: %s",path.c_str(),strerror(errno));
//		exit(0);
//	}
//	assert(mt_fs);
//
//	size_t dsize = 0;
//	size_t index = 0;
//	while(index<capacity&&fread((void *)&dsize, sizeof(size_t), 1, mt_fs)>0){
//		HiMesh_Wrapper *w = new HiMesh_Wrapper();
//		w->offset = dsize;
//		fread((void *)&w->data_size, sizeof(size_t), 1, mt_fs);
//		w->id = index++;
//		w->box.id = w->id;
//		// read the voxels into the wrapper
//		fread((void *)&dsize, sizeof(size_t), 1, mt_fs);
//		for(int i=0;i<dsize;i++){
//			Voxel *v = new Voxel();
//			fread((void *)v->low, sizeof(float), 3, mt_fs);
//			fread((void *)v->high, sizeof(float), 3, mt_fs);
//			fread((void *)v->core, sizeof(float), 3, mt_fs);
//			w->voxels.push_back(v);
//			w->box.update(*v);
//		}
//		objects.push_back(w);
//		space.update(w->box);
//	}
//	fclose(mt_fs);
//	return true;
//}

}


