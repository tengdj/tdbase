/*
 * tile_IO.cpp
 *
 *  Created on: Apr 4, 2023
 *      Author: teng
 */

#include "tile.h"

namespace tdbase{

// do the initialization job
void Tile::load(){
	struct timeval start = get_cur_time();
	if(!file_exist(tile_path.c_str())){
		log("%s does not exist", tile_path.c_str());
		exit(-1);
	}
	// load the raw data into the buffer
	data_size = file_size(tile_path.c_str());
	data_buffer = new char[data_size];
	//process_lock();
	FILE *dt_fs = fopen(tile_path.c_str(), "r");
	if(dt_fs == NULL){
		log("failed to open file %s",tile_path.c_str());
	}
	if(fread((void *)data_buffer, sizeof(char), data_size, dt_fs) != data_size){
		assert("failed reading" && false);
	}
	fclose(dt_fs);
	//process_unlock();

	// parsing the metadata from the dt file
	Decoding_Type dtype = (Decoding_Type)data_buffer[0];
	size_t offset = 1;// the first byte is the file type, raw or compressed
	size_t index = 0;
	while(offset < data_size){
		// create a wrapper with the meta information
		HiMesh_Wrapper * w = new HiMesh_Wrapper(data_buffer + offset, index++, dtype);
		offset += w->data_size + w->meta_size + sizeof(size_t);
		objects.push_back(w);
		space.update(w->box);
	}

	tree = build_octree(10);
	logt("loaded %ld polyhedra in tile %s", start, objects.size(), tile_path.c_str());
}

void Tile::dump_compressed(const char *path){
	ofstream *os = new std::ofstream(path, std::ios::out | std::ios::binary);
	assert(os);
	char type = (char)COMPRESSED;
	os->write(&type, 1);
	for(HiMesh_Wrapper *wr:objects){
		assert(wr->type == COMPRESSED);
		HiMesh *nmesh = wr->get_mesh();
		//tdbase::write_polyhedron(&shifted, ids++);
		size_t size = nmesh->get_data_size();
		os->write((char *)&size, sizeof(size_t));
		os->write(nmesh->get_data(), nmesh->get_data_size());
		size = wr->voxels.size();
		os->write((char *)&size, sizeof(size_t));
		for(Voxel *v:wr->voxels){
			os->write((char *)v->low, 3*sizeof(float));
			os->write((char *)v->high, 3*sizeof(float));
			os->write((char *)v->core, 3*sizeof(float));
		}
	}
	os->close();
}

// dump to a raw format tile file
void Tile::dump_raw(const char *path){

	ofstream *os = new std::ofstream(path, std::ios::out | std::ios::binary);

	char *buffer = new char[data_size*20];
	size_t offset = 0;
	buffer[0] = (char)RAW;
	offset++;

	for(HiMesh_Wrapper *wr:objects){
		assert(wr->type != RAW && "already be in raw format");

		size_t *dsize_holder = (size_t *)(buffer+offset);
		offset += sizeof(size_t);

		const size_t st_offset = offset;

		map<int, float> hausdorffs;
		map<int, float> proxyhausdorffs;
		for(int lod=20;lod<=100;lod+=20){
			wr->decode_to(lod);
				if(lod!=100){
				hausdorffs[lod] = wr->get_mesh()->collectGlobalHausdorff().second;
				proxyhausdorffs[lod] = wr->get_mesh()->collectGlobalHausdorff().first;
			}else{
				hausdorffs[lod] = 0;
				proxyhausdorffs[lod] = 0;
			}

			for(Voxel *v:wr->voxels){
				v->offset_lod[lod] = offset - st_offset;
				v->volume_lod[lod] = v->num_triangles;
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
				*(size_t *)(buffer+offset) = v->volume_lod[lod];
				offset += sizeof(size_t);
			}
		}
	}

	os->write(buffer, offset);
	os->close();
	delete os;
	delete []buffer;
	log("converted to %s",path);
}

void Tile::dump_sql(const char *path, const char *table){
	std::filebuf fb;
	fb.open(path, std::ios::out | std::ios::trunc);
	if(fb.is_open()) {
		std::ostream os(&fb);
		for(HiMesh_Wrapper *wr:objects){
			string wkt = wr->get_mesh()->to_wkt();
			os << "INSERT INTO "<<table<<"(id, hausdorff, phausdorff, geom) VALUES ("<<wr->id<<","<<wr->getHausdorffDistance()<<","<<wr->getProxyHausdorffDistance()<<",ST_GeomFromEWKT('"<<wkt.c_str()<<"'));"<<endl;
		}
		fb.close();
	}else{
		std::cerr<<"cannot find path "<<path<<std::endl;
	}
}

}


