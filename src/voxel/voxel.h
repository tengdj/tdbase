/*
 * voxel.h
 *
 *  Created on: Nov 11, 2019
 *      Author: teng
 */

#ifndef VOXEL_H_
#define VOXEL_H_
#include <vector>
#include <queue>

#include "../PPMC/mymesh.h"


namespace hispeed{

enum voxel_type{
	VT_EDGE,
	VT_TRIANGLE
};


class voxel{
	// id
	long id = 0;
	// mbb
	mbb box;
	enum voxel_type type = VT_EDGE;
	int data_length[10];
	// the decoded edges, or triangles
	char *data[10];

public:

	voxel(){
		for(int i=0;i<10;i++){
			data[i] = NULL;
		}
	};

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


}




#endif /* VOXEL_H_ */
