/*
 * tile.h
 *
 *  Created on: Nov 15, 2019
 *      Author: teng
 */

#ifndef HISPEED_TILE_H_
#define HISPEED_TILE_H_

namespace hispeed{

enum voxel_type{
	VT_EDGE,
	VT_TRIANGLE
};

class voxel{
	// id
	long id = 0;
	// boundary box of the voxel
	aab box;
	enum voxel_type type = VT_EDGE;
	int data_length[10];
	// the decoded edges, or triangles
	char *data[10];
	MyMesh *mesh = NULL;

public:

	voxel(){
		for(int i=0;i<10;i++){
			data[i] = NULL;
		}
	};
	~voxel(){
		for(int i=0;i<10;i++){
			if(data[i]!=NULL){
				delete data[i];
				data[i] = NULL;
			}
		}
	}

	bool decompressed(int lod){
		assert(lod>=0&&lod<10);
		return data[lod]!=NULL;
	}

	float *get_edges(int lod){
		assert(lod>=0&&lod<10);
		return (float *)data[lod];
	}

	float *get_triangles(int lod){
		assert(lod>=0&&lod<10);
		return (float *)data[lod];
	}


};

typedef struct voxel_pair{
	voxel *A;
	voxel *B;
} voxel_pair;


class tile{

	std::vector<voxel *> voxels;
public:

	~tile(){
		for(voxel *v:voxels){
			delete v;
		}
	}


};


}



#endif /* HISPEED_TILE_H_ */
