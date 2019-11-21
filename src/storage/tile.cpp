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
// tile->voxel_groups->voxels
Tile::Tile(std::string path){
	set_prefix(path);
	active = load();
	dt_fs = fopen(data_path.c_str(), "r");
	active &= (dt_fs!=NULL);
	// has to be true
	assert(active);
}

Tile::~Tile(){
	for(Voxel_group *v:voxel_groups){
		delete v;
	}
	// close the data file pointer
	if(dt_fs!=NULL){
		fclose(dt_fs);
	}
}

void Tile::set_prefix(string path){
	meta_path = path+".mt";
	data_path = path+".dt";
}

void Tile::add_polyhedron(MyMesh *mesh, long offset){
	assert(mesh && mesh->i_mode == DECOMPRESSION_MODE_ID);
	Voxel_group *g = new Voxel_group();
	g->offset = offset;
	g->data_size = mesh->dataOffset;
	mesh->generate_mbbs();
	vector<aab> mbbs = mesh->get_mbbs();
	for(aab b:mbbs){
		g->add_voxel(b);
	}
	add_group(g);
}

// parse the meta data from the compressed mesh data
bool Tile::parse_raw(){
	cerr<<"parsing from "<<data_path<<endl;
	assert(!active);
	FILE* fs = fopen(data_path.c_str(), "r");
	if(fs==NULL){
		cerr<<"failed to open file "<<meta_path<<endl;
		return false;
	}
	long fsize = hispeed::file_size(data_path.c_str());
	// must be empty
	assert(voxel_groups.size()==0);
	size_t dsize = 0;
	char *buf = new char[3*1024*1024];
	aab tmp;
	long offset = 0;
	int next_report = 0;
	int num_objects = 0;
	while(fread((void *)&dsize, sizeof(size_t), 1, fs)>0){
		assert(dsize<3*1024*1024);
		// read data
		assert(dsize==fread((void *)buf, sizeof(char), dsize, fs));
		// decompress the data and add it to the Tile
		MyMesh *mesh = hispeed::compress_mesh(buf, dsize, true);
		add_polyhedron(mesh, offset);
		delete mesh;
		// update the offset for next
		offset += dsize+sizeof(size_t);
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
 * to disk
 *
 * */
bool Tile::persist(){
	FILE* fs = fopen(meta_path.c_str(), "wb+");
	if(fs==NULL){
		cerr<<"failed to open file "<<meta_path<<endl;
		return false;
	}

	fwrite(&id, 1, sizeof(int), fs);
	uint size = voxel_groups.size();
	fwrite(&size, 1, sizeof(uint), fs);
	for(Voxel_group *g:voxel_groups){
		if(!g->persist(fs)){
			return false;;
		}
	}
	fclose(fs);
	active = true;
	return true;
}

/*
 * load the persisted information from disk
 * to construct the tile. Note that the information
 * of the voxels in the tile are parsed lazily only
 * on-demand. Thus the load function only load the description
 * info. the compressed data will be loaded separately
 * from another file during query.
 * */
bool Tile::load(){
	FILE* fs = fopen(meta_path.c_str(), "r");
	if(fs==NULL){
		cerr<<"failed to open file "<<meta_path<<endl;
		return false;
	}

	fread((void *)&id, sizeof(int), 1, fs);
	uint size = 0;
	fread((void *)&size, sizeof(uint), 1, fs);
	for(uint i=0;i<size;i++){
		Voxel_group *g = new Voxel_group();
		if(!g->load(fs)){
			return false;;
		}
		add_group(g);
	}
	fclose(fs);
	active = true;
	return true;
}


// retrieve the mesh of the voxel group with ID id on demand
void Tile::retrieve_mesh(int id){
	assert(id>=0&&id<voxel_groups.size());
	Voxel_group *cur_group = voxel_groups[id];
	assert(cur_group && cur_group->mesh==NULL);
	cur_group->mesh_data = new char[cur_group->data_size];
	assert(dt_fs);
	assert(cur_group->data_size==
			fread(cur_group->mesh_data, sizeof(char), cur_group->data_size, dt_fs));
	cur_group->mesh =
			hispeed::compress_mesh(cur_group->mesh_data, cur_group->data_size);
}

// release the space for mesh object and the triangles
// decoded in each voxel groups
bool Tile::release_mesh(int id){
	assert(id>=0&&id<this->voxel_groups.size());
	voxel_groups[id]->release_data();
	return true;
}


void Tile::print(){
	cout<<"tile "<<id<<endl;
	for(Voxel_group *g:voxel_groups){
		g->print();
	}
	cout<<endl;
}

}

