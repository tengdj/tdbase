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
Tile::Tile(std::string path, size_t capacity){
	struct timeval start = get_cur_time();
	this->capacity = capacity;
	if(!hispeed::file_exist(path.c_str())){
		log("%s does not exist", path.c_str());
		exit(-1);
	}
	dt_fs = fopen(path.c_str(), "r");
	if(!dt_fs){
		log("%s can not be opened", path.c_str());
		exit(-1);
	}
	string meta_path = path;
	boost::replace_all(meta_path, ".dt", ".mt");
	if(!hispeed::file_exist(meta_path.c_str())){
		parse_raw();
		persist(meta_path);
	}else{
		load(meta_path);
	}
	pthread_mutex_init(&lock, NULL);
	logt("loaded %ld polyhedra in tile %s", start, objects.size(), path.c_str());
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

// persist the meta data for current tile as a cache
bool Tile::persist(string path){
	FILE *mt_fs = fopen(path.c_str(), "wb+");
	assert(mt_fs);
	for(HiMesh_Wrapper *w:objects){
		fwrite((void *)&w->offset, sizeof(size_t), 1, mt_fs);
		fwrite((void *)&w->data_size, sizeof(size_t), 1, mt_fs);
		size_t dsize = w->voxels.size();
		fwrite((void *)&dsize, sizeof(size_t), 1, mt_fs);
		for(Voxel *v:w->voxels){
			fwrite((void *)v->box.min, sizeof(float), 3, mt_fs);
			fwrite((void *)v->box.max, sizeof(float), 3, mt_fs);
			fwrite((void *)v->core, sizeof(float), 3, mt_fs);
		}
	}
	fclose(mt_fs);
	return true;
}

// load from the cached data
bool Tile::load(string path){
	FILE *mt_fs = fopen(path.c_str(), "r");
	assert(mt_fs);
	size_t dsize = 0;
	size_t index = 0;
	while(index<capacity&&fread((void *)&dsize, sizeof(size_t), 1, mt_fs)>0){
		HiMesh_Wrapper *w = new HiMesh_Wrapper();
		w->offset = dsize;
		fread((void *)&w->data_size, sizeof(size_t), 1, mt_fs);
		w->id = index++;
		w->box.id = w->id;
		// read the voxels into the wrapper
		fread((void *)&dsize, sizeof(size_t), 1, mt_fs);
		for(int i=0;i<dsize;i++){
			Voxel *v = new Voxel();
			fread((void *)v->box.min, sizeof(float), 3, mt_fs);
			fread((void *)v->box.max, sizeof(float), 3, mt_fs);
			fread((void *)v->core, sizeof(float), 3, mt_fs);
			w->voxels.push_back(v);
			w->box.box.update(v->box);
		}
		objects.push_back(w);
		box.update(w->box.box);
	}
	fclose(mt_fs);
	return true;
}

// parse the metadata
bool Tile::parse_raw(){
	assert(dt_fs);
	size_t dsize = 0;
	long offset = 0;
	size_t index = 0;
	fseek(dt_fs, 0, SEEK_SET);
	while(index<capacity&&fread((void *)&dsize, sizeof(size_t), 1, dt_fs)>0){
		fseek(dt_fs, dsize, SEEK_CUR);
		HiMesh_Wrapper *w = new HiMesh_Wrapper();
		offset += sizeof(size_t);
		w->offset = offset;
		w->data_size = dsize;
		w->id = index++;
		w->box.id = w->id;
		// read the voxels into the wrapper
		fread((void *)&dsize, sizeof(size_t), 1, dt_fs);
		for(int i=0;i<dsize;i++){
			Voxel *v = new Voxel();
			fread((void *)v->box.min, sizeof(float), 3, dt_fs);
			fread((void *)v->box.max, sizeof(float), 3, dt_fs);
			fread((void *)v->core, sizeof(float), 3, dt_fs);
			w->voxels.push_back(v);
			w->box.box.update(v->box);
		}
		objects.push_back(w);
		box.update(w->box.box);
		// update the offset for next
		offset += w->data_size+sizeof(size_t)+9*sizeof(float)*dsize;
	}
	return true;
}

// retrieve the mesh of the voxel group with ID id on demand
void Tile::retrieve_mesh(int id){
	assert(id>=0&&id<objects.size());
	HiMesh_Wrapper *wrapper = objects[id];
	if(wrapper->mesh==NULL){
		char *mesh_data = new char[wrapper->data_size];
		assert(dt_fs);
		fseek(dt_fs, wrapper->offset, SEEK_SET);
		assert(wrapper->data_size==
				fread(mesh_data, sizeof(char), wrapper->data_size, dt_fs));
		wrapper->mesh = new HiMesh(mesh_data, wrapper->data_size);
		delete mesh_data;
	}
}

OctreeNode *Tile::build_octree(size_t leaf_size){
	OctreeNode *octree = new OctreeNode(box, 0, leaf_size);
	for(HiMesh_Wrapper *w:objects){
		octree->addObject(&w->box);
	}
	return octree;
}

void Tile::decode_to(int id, int lod){
	assert(id>=0&&id<objects.size());
	pthread_mutex_lock(&lock);
	if(objects[id]->mesh==NULL){
		retrieve_mesh(id);
	}
	assert(objects[id]->mesh);
	objects[id]->mesh->advance_to(lod);
	pthread_mutex_unlock(&lock);
}

}

