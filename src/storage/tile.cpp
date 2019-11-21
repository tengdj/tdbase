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
	set_path(path);
	active = load()&&
			(dt_fs = fopen(data_path.c_str(), "r"));
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

void Tile::set_path(string path){
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

/*
 * persist the information of the tile
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
	fwrite(&space.min, 3, sizeof(float), fs);
	fwrite(&space.max, 3, sizeof(float), fs);
	uint size = voxel_groups.size();
	fwrite(&size, 1, sizeof(uint), fs);
	for(Voxel_group *g:voxel_groups){
		if(!g->persist(fs)){
			return false;;
		}
	}
	fclose(fs);
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

	fread((void *)&id, 1, sizeof(int), fs);
	fread((void *)&space.min, 3, sizeof(float), fs);
	fread((void *)&space.max, 3, sizeof(float), fs);
	uint size = 0;
	fread((void *)&size, 1, sizeof(uint), fs);
	for(uint i=0;i<size;i++){
		Voxel_group *v = new Voxel_group();
		v->id = i;
		v->tile = this;
		if(!v->load(fs)){
			return false;;
		}
	}
	fclose(fs);
	return true;
}


// retrieve the mesh of the voxel group with ID id on demand
void Tile::retrieve_mesh(int id){
	assert(id>=0&&id<voxel_groups.size());
	Voxel_group *cur_group = voxel_groups[id];
	assert(cur_group && cur_group->mesh==NULL);
	cur_group->mesh_data = new char[cur_group->data_size];
	assert(dt_fs);
	fread(cur_group->mesh_data, cur_group->data_size, sizeof(char), dt_fs);
	cur_group->mesh =
			hispeed::get_compressed_mesh(cur_group->mesh_data, cur_group->data_size);
}

// release the space for mesh object and the triangles
// decoded in each voxel groups
bool Tile::release_mesh(int id){
	assert(id>=0&&id<this->voxel_groups.size());
	voxel_groups[id]->release_data();
	return true;
}

}

