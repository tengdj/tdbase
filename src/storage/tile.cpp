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

namespace tdbase{

// load meta data from file and construct the hierarchy structure
Tile::Tile(std::string path, size_t capacity, bool active_load){
	tile_path = path;
	tile_capacity = capacity;
	if(active_load){
		load();
	}
}

Tile::Tile(std::vector<HiMesh_Wrapper *> &objs){
	objects.assign(objs.begin(), objs.end());
	for(HiMesh_Wrapper *wr:objs){
		if(wr->type == MULTIMESH){
			for(auto mesh:wr->get_meshes()){
				data_size += sizeof(float)*3*mesh.second->size_of_triangles();
			}
		}
	}

#pragma omp parallel for
	for(HiMesh_Wrapper *wr:objs){
		if(wr->type == MULTIMESH){
			HiMesh *original = wr->get_meshes()[100];
			for(int lod=20;lod<=80;lod+=20){
				HiMesh *low = wr->get_meshes()[lod];
				low->computeHausdorfDistance(original);
				//pair<float, float> hf = low->collectGlobalHausdorff(AVG);
				//log("%d: %f %f", lod, hf.first, hf.second);
			}
			//log("");
		}
	}
}

Tile::~Tile(){
	for(HiMesh_Wrapper *h:objects){
		delete h;
	}
	if(data_buffer!=NULL){
		delete []data_buffer;
	}
	if(tree){
		delete tree;
	}
}

HiMesh *Tile::get_mesh(int id){
	return objects[id]->get_mesh();
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

