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
Tile::Tile(std::string prefix){
	this->prefix = prefix;
	if(!hispeed::file_exist(get_data_path().c_str())){
		cerr<<get_data_path()<<" does not exist"<<endl;
		exit(-1);
	}
	if(!hispeed::file_exist(get_meta_path().c_str())){
		parse_raw();
		persist();
	}else{
		load();
	}
	dt_fs = fopen(get_data_path().c_str(), "r");
}

Tile::~Tile(){
	for(HiMesh_Wrapper *h:objects){
		delete h;
	}
	// close the data file pointer
	if(dt_fs!=NULL){
		fclose(dt_fs);
	}
}



// parse the meta data from the compressed mesh data
bool Tile::parse_raw(){
	string data_path = get_data_path();
	cerr<<"parsing from "<<data_path<<endl;
	FILE* fs = fopen(data_path.c_str(), "r");
	if(fs==NULL){
		cerr<<"failed to open file "<<data_path<<endl;
		return false;
	}
	long fsize = hispeed::file_size(data_path.c_str());
	// must be empty
	assert(objects.size()==0);
	size_t dsize = 0;
	char *buf = new char[3*1024*1024];
	aab tmp;
	long offset = 0;
	int next_report = 0;
	int num_objects = 0;
	int index = 0;
	while(fread((void *)&dsize, sizeof(size_t), 1, fs)>0){
		assert(dsize<3*1024*1024);
		// read data
		assert(dsize==fread((void *)buf, sizeof(char), dsize, fs));
		// decompress the data and add it to the Tile
		MyMesh *mesh = hispeed::decompress_mesh(buf, dsize, true);
		assert(dsize == mesh->dataOffset);
		HiMesh_Wrapper *w = new HiMesh_Wrapper();
		offset += sizeof(size_t);
		w->offset = offset;
		w->data_size = dsize;
		w->id = index++;
		w->box = aab(mesh->bbMin[0], mesh->bbMin[1], mesh->bbMin[2],
					 mesh->bbMax[0], mesh->bbMax[1], mesh->bbMax[2]);
		objects.push_back(w);
		delete mesh;
		// update the offset for next
		offset += dsize;
		num_objects++;
		if(offset*100/fsize==next_report){
			cerr<<"processed "<<num_objects<<" objects\t("<<next_report<<"%)"<<endl;
			next_report++;
		}
	}
	return true;
}

/*
 * persist the meta data of the tile
 * to disk.
 *
 * */
bool Tile::persist(){

	FILE* fs = fopen(get_meta_path().c_str(), "wb+");
	if(fs==NULL){
		cerr<<"failed to open file "<<get_meta_path()<<endl;
		return false;
	}
	size_t size = objects.size();
	fwrite(&size, sizeof(size_t), 1, fs);
	for(HiMesh_Wrapper *w:objects){
		fwrite((void *)w->box.min, sizeof(float), 3, fs);
		fwrite((void *)w->box.max, sizeof(float), 3, fs);
		fwrite((void *)&w->offset, sizeof(long), 1, fs);
		fwrite((void *)&w->data_size, sizeof(long), 1, fs);
	}
	fclose(fs);
	return true;
}

/*
 * load the persisted information from disk
 * to construct the tile. The load function only load the description
 * info of each mesh. the compressed data will be loaded separately
 * from another file during query.
 * */
bool Tile::load(){
	string meta_path = get_meta_path();
	FILE* fs = fopen(meta_path.c_str(), "r");
	if(fs==NULL){
		cerr<<"failed to open file "<<meta_path<<endl;
		return false;
	}
	size_t size = 0;
	fread((void *)&size, sizeof(size_t), 1, fs);
	for(int i=0;i<size;i++){
		HiMesh_Wrapper *h = new HiMesh_Wrapper();
		h->id = i;
		fread((void *)h->box.min, sizeof(float), 3, fs);
		fread((void *)h->box.max, sizeof(float), 3, fs);
		fread((void *)&h->offset, sizeof(long), 1, fs);
		fread((void *)&h->data_size, sizeof(long), 1, fs);
		objects.push_back(h);
	}
	fclose(fs);
	return true;
}


// retrieve the mesh of the voxel group with ID id on demand
void Tile::retrieve_mesh(int id){
	assert(id>=0&&id<objects.size());
	HiMesh_Wrapper *wrapper = objects[id];
	assert(wrapper->mesh==NULL);
	char *mesh_data = new char[wrapper->data_size];
	assert(dt_fs);
	fseek(dt_fs, wrapper->offset, SEEK_SET);
	assert(wrapper->data_size==
			fread(mesh_data, sizeof(char), wrapper->data_size, dt_fs));
	wrapper->mesh = new HiMesh(mesh_data, wrapper->data_size);
	delete mesh_data;
}


void Tile::print(){
	int index = 0;
	for(HiMesh_Wrapper *w:objects){
		cout<<index++<<"\t"<<w->box<<endl;
	}
	cout<<endl;
}

}

