/*
 * joiner.cpp
 *
 *  Created on: Nov 20, 2019
 *      Author: teng
 */




#include "../storage/tile.h"
#include "../geometry/geometry.h"
using namespace std;
using namespace hispeed;

int main(int argc, char **argv){
//
//	int lod = 100;
//	if(argc>1){
//		lod = atoi(argv[1]);
//	}
//
//	int totalnum = 10000;
//	if(argc>2){
//		totalnum = atoi(argv[2]);
//	}
//	struct timeval start = get_cur_time();
//	Tile *tile = new Tile("test");
//	printf("loading tile meta takes %f milliseconds\n", get_time_elapsed(start, true));
//	for(int i=0;i<tile->num_objects();i++){
//		tile->get_mesh(i, lod);
//	}
//	printf("decompressing takes %f milliseconds\n", get_time_elapsed(start, true));
//
//	HiMesh *geom1 = tile->get_mesh(0, lod);
//	HiMesh *geom2 = tile->get_mesh(1, lod);
//
//	size_t size1 = geom1->get_segment_num();
//	size_t size2 = geom2->get_segment_num();
//	float *data1 = geom1->get_segments();
//	float *data2 = geom2->get_segments();
//	printf("%ld X %ld = %ld \n", size1, size2, size1*size2);
//	printf("loading segments takes %f milliseconds\n", get_time_elapsed(start, true));
//
//	size_t voxel_size = 400;
//	voxel_size = std::min(voxel_size, size1);
//	voxel_size = std::min(voxel_size, size2);
//
//	float *wdata1 = new float[6*voxel_size*totalnum];
//	float *wdata2 = new float[6*voxel_size*totalnum];
//	for(int i=0;i<totalnum;i++){
//		memcpy((char *)(wdata1+i*6*voxel_size), (char *)data1, voxel_size*6*sizeof(float));
//		memcpy((char *)(wdata2+i*6*voxel_size), (char *)data2, voxel_size*6*sizeof(float));
//	}
//	printf("copy data takes %f milliseconds\n", get_time_elapsed(start, true));
//	float *distances = new float[totalnum];
//	for(int i=0;i<totalnum;i++){
//		distances[i] = SegDist_single(wdata1+voxel_size*6*i, wdata2+voxel_size*6*i, voxel_size, voxel_size);
//	}
//	printf("get distance %f with CPU take %f miliseconds\n",std::sqrt(distances[0]),get_time_elapsed(start, true));
//
//	hispeed::SegDist_batch_gpu(wdata1, wdata2, voxel_size, totalnum, distances);
//	printf("get distance %f with GPU take %f miliseconds\n",std::sqrt(distances[0]),get_time_elapsed(start, true));
//
//	delete tile;
//	delete data1;
//	delete data2;
//	delete wdata1;
//	delete wdata2;
//	delete distances;
}
