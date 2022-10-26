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
Tile::Tile(std::string path, size_t capacity, bool init_data){
	if(!hispeed::file_exist(path.c_str())){
		log("%s does not exist", path.c_str());
		exit(-1);
	}
	tile_path = path;
	tile_capacity = capacity;
	if(init_data){
		init();
	}
}

void Tile::init(){
	struct timeval start = get_cur_time();
	FILE *dt_fs = fopen(tile_path.c_str(), "r");
	if(!dt_fs){
		log("%s can not be opened", tile_path.c_str());
		exit(-1);
	}
	string meta_path = tile_path;
	boost::replace_all(meta_path, ".dt", ".mt");
	process_lock();
	if(!hispeed::file_exist(meta_path.c_str())){
		parse_raw(dt_fs);
		persist(meta_path);
	}else{
		load(meta_path, tile_capacity);
	}
	size_t sz = file_size(tile_path.c_str());
	data_buffer = new char[sz];
	fseek(dt_fs, 0, SEEK_SET);
	size_t rdsze = fread((void *)data_buffer, sizeof(char), sz, dt_fs);
	assert(rdsze == sz);
	process_unlock();

//	size_t sm = 0;
//	for(int i=0;i<sz;i++){
//		sm += data_buffer[i];
//	}
//	log("%ld",sm);
	// close the data file pointer
	fclose(dt_fs);
	if(!global_ctx.use_multimbb){
		disable_innerpart();
	}
	logt("loaded %ld polyhedra in tile %s", start, objects.size(), tile_path.c_str());
}

Tile::~Tile(){
	for(HiMesh_Wrapper *h:objects){
		delete h;
	}
	if(data_buffer!=NULL){
		delete []data_buffer;
	}
}

void Tile::disable_innerpart(){
	for(HiMesh_Wrapper *w:this->objects){
		if(w->voxels.size()>1){
			for(Voxel *v:w->voxels){
				delete v;
			}
			w->voxels.clear();
			Voxel *v = new Voxel();
			v->set_box(w->box);
			for(int i=0;i<3;i++){
				v->core[i] = (v->high[i]-v->low[i])/2+v->low[i];
			}
			w->voxels.push_back(v);
		}
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
			fwrite((void *)v->low, sizeof(float), 3, mt_fs);
			fwrite((void *)v->high, sizeof(float), 3, mt_fs);
			fwrite((void *)v->core, sizeof(float), 3, mt_fs);
		}
	}
	fclose(mt_fs);
	return true;
}

// load from the cached data
bool Tile::load(string path, int capacity){
	FILE *mt_fs = fopen(path.c_str(), "r");
	if(mt_fs==NULL){
		log("%s cannot be opened, error: %s",path.c_str(),strerror(errno));
		exit(0);
	}
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
			fread((void *)v->low, sizeof(float), 3, mt_fs);
			fread((void *)v->high, sizeof(float), 3, mt_fs);
			fread((void *)v->core, sizeof(float), 3, mt_fs);
			w->voxels.push_back(v);
			w->box.update(*v);
		}
		objects.push_back(w);
		box.update(w->box);
	}
	fclose(mt_fs);
	return true;
}

// parse the metadata
bool Tile::parse_raw(FILE *dt_fs){
	assert(dt_fs);
	size_t dsize = 0;
	long offset = 0;
	size_t index = 0;
	fseek(dt_fs, 0, SEEK_SET);
	while(fread((void *)&dsize, sizeof(size_t), 1, dt_fs)>0){
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
			fread((void *)v->low, sizeof(float), 3, dt_fs);
			fread((void *)v->high, sizeof(float), 3, dt_fs);
			fread((void *)v->core, sizeof(float), 3, dt_fs);
			w->voxels.push_back(v);
			w->box.update(*v);
		}
		objects.push_back(w);
		box.update(w->box);
		// update the offset for next
		offset += w->data_size+sizeof(size_t)+9*sizeof(float)*dsize;
	}
	return true;
}

// retrieve the mesh of the voxel group with ID id on demand
void Tile::retrieve_mesh(size_t id){

	assert(id>=0&&id<objects.size());
	HiMesh_Wrapper *wrapper = objects[id];
	if(wrapper->mesh==NULL){
		timeval cur = hispeed::get_cur_time();
		wrapper->mesh = new HiMesh(data_buffer+wrapper->offset, wrapper->data_size);
		newmesh_time += hispeed::get_time_elapsed(cur, true);
	}
	assert(wrapper->mesh);
}

OctreeNode *Tile::build_octree(size_t leaf_size){
	OctreeNode *octree = new OctreeNode(box, 0, leaf_size);
	for(HiMesh_Wrapper *w:objects){
		octree->addObject(&w->box);
	}
	return octree;
}


void Tile::decode_to(size_t id, uint lod){
	assert(lod>=0 && lod<=100);
	assert(id>=0&&id<objects.size());
	timeval cur = hispeed::get_cur_time();
	timeval start = hispeed::get_cur_time();
	retrieve_mesh(id);
	retrieve_time += hispeed::get_time_elapsed(cur,true);
	objects[id]->advance_to(lod);
	advance_time += hispeed::get_time_elapsed(cur,true);
	decode_time += hispeed::get_time_elapsed(start,true);
}

void Tile::add_raw(char *data){
	size_t offset = 0;
	size_t size_tmp = 0;
	memcpy((char *)&size_tmp, data+offset, sizeof(size_t));
	offset += sizeof(size_t);
	HiMesh *mesh = new HiMesh(data+offset, size_tmp);
	offset += size_tmp;
	memcpy((char *)&size_tmp, data+offset, sizeof(size_t));
	offset += sizeof(size_t);
	HiMesh_Wrapper *hw = new HiMesh_Wrapper();
	for(size_t i=0;i<size_tmp;i++){
		Voxel *v = new Voxel();
		memcpy((char *)v->low, data+offset, 3*sizeof(float));
		offset += 3*sizeof(float);
		memcpy((char *)v->high, data+offset, 3*sizeof(float));
		offset += 3*sizeof(float);
		memcpy((char *)v->core, data+offset, 3*sizeof(float));
		offset += 3*sizeof(float);
		hw->voxels.push_back(v);
		hw->box.update(*v);
	}
	this->box.update(hw->box);

}

}

