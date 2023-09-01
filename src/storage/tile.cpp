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

// load meta data from file and construct the hierarchy structure
Tile::Tile(std::string path, size_t capacity, Decoding_Type dt){
	dtype = dt;
	tile_path = path;
	tile_capacity = capacity;
	init();
}

Tile::~Tile(){
	for(HiMesh_Wrapper *h:objects){
		delete h;
	}
	if(data_buffer!=NULL){
		delete []data_buffer;
	}
}

// do the initialization job
void Tile::init(){
	struct timeval start = get_cur_time();
	if(!file_exist(tile_path.c_str())){
		log("%s does not exist", tile_path.c_str());
		exit(-1);
	}

	process_lock();
	// load the raw data into the buffer
	data_size = file_size(tile_path.c_str());
	data_buffer = new char[data_size];
	FILE *dt_fs = fopen(tile_path.c_str(), "r");
	assert(fread((void *)data_buffer, sizeof(char), data_size, dt_fs) == data_size);
	fclose(dt_fs);
	process_unlock();

	// parsing the metadata from the dt file
	size_t offset = 0;
	size_t index = 0;
	while(offset < data_size){

		// create a wrapper with the meta information
		HiMesh_Wrapper * w = new HiMesh_Wrapper(data_buffer + offset, index++, dtype);
		offset += w->data_size + w->meta_size + sizeof(size_t);
		objects.push_back(w);
		space.update(w->box);
	}

	logt("loaded %ld polyhedra in tile %s", start, objects.size(), tile_path.c_str());
}

void Tile::convert_raw(const char *path){
	assert(dtype == COMPRESSED);
	size_t offset = 0;
	char *buffer = new char[data_size*20];

	for(HiMesh_Wrapper *wr:objects){
		size_t *dsize_holder = (size_t *)(buffer+offset);
		offset += sizeof(size_t);

		const size_t st_offset = offset;

		map<int, float> hausdorffs;
		map<int, float> proxyhausdorffs;
		for(int lod=20;lod<=100;lod+=20){
			wr->decode_to(lod);
			hausdorffs[lod] = wr->mesh->getHausdorffDistance();
			proxyhausdorffs[lod] = wr->mesh->getProxyHausdorffDistance();

			for(Voxel *v:wr->voxels){
				v->offset_lod[lod] = offset - st_offset;
				v->volumn_lod[lod] = v->num_triangles;
				memcpy(buffer+offset, v->triangles, v->num_triangles*sizeof(float)*9);
				offset += v->num_triangles*sizeof(float)*9;
				memcpy(buffer+offset, v->hausdorff, v->num_triangles*sizeof(float)*2);
				offset += v->num_triangles*sizeof(float)*2;
			}
		}
		// update the data size
		*dsize_holder = offset-st_offset;

		// store the voxel number (put it here for aligning with the decoding mode)
		*(size_t *)(buffer + offset) = wr->voxels.size();
		offset += sizeof(size_t);

		// store the polyhedron-level hausdorff information for all the LODs
		for(int lod=20;lod<=100;lod+=20){
			*(float *)(buffer + offset) = hausdorffs[lod];
			offset += sizeof(float);
			*(float *)(buffer + offset) = proxyhausdorffs[lod];
			offset += sizeof(float);
		}

		// store the voxel information
		for(Voxel *v:wr->voxels){
			memcpy(buffer+offset, v->low, sizeof(float)*3);
			offset += 3*sizeof(float);
			memcpy(buffer+offset, v->high, sizeof(float)*3);
			offset += 3*sizeof(float);
			memcpy(buffer+offset, v->core, sizeof(float)*3);
			offset += 3*sizeof(float);
			for(int lod=20;lod<=100;lod+=20){
				*(size_t *)(buffer+offset) = v->offset_lod[lod];
				offset += sizeof(size_t);
				*(size_t *)(buffer+offset) = v->volumn_lod[lod];
				offset += sizeof(size_t);
			}
		}
	}

	ofstream *os = new std::ofstream(path, std::ios::out | std::ios::binary);
	os->write(buffer, offset);
	os->close();
	delete os;
}

HiMesh *Tile::get_mesh(int id){
	return objects[id]->mesh;
}

void Tile::decode_all(int lod){
	for(HiMesh_Wrapper *w:objects){
		w->decode_to(lod);
	}
}

OctreeNode *Tile::build_octree(size_t leaf_size){
	OctreeNode *octree = new OctreeNode(space, 0, leaf_size);
	for(HiMesh_Wrapper *w:objects){
		octree->addObject(&w->box);
	}
	return octree;
}

}

