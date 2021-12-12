/*
 * decomp.cpp
 *
 *  Created on: Oct 22, 2019
 *      Author: teng
 */


#include "../PPMC/ppmc.h"
#include "../util/util.h"
#include "../spatial/spatial.h"
#include "../spatial/himesh.h"
#include <algorithm>
#include "zlib.h"

using namespace hispeed;

int main(int argc, char **argv){


	int start_lod = 0;
	int end_lod = 10;
	if(argc>1){
		start_lod = atoi(argv[1]);
		end_lod = start_lod;
	}
	assert(start_lod>=0&&start_lod<=10);
	// Init the random number generator.
	log("start compressing");
	string mesh_str = read_off_stdin();
	MyMesh *compressed = get_mesh(mesh_str);

	struct timeval starttime = get_cur_time();
	//assert(compressed->size_of_border_edges()&&"must be manifold");
	compressed->completeOperation();
	logt("compress", starttime);

	MyMesh *testc[100];
	for(int i=0;i<100;i++){
		testc[i] = get_mesh(mesh_str);
	}
	struct timeval sst = get_cur_time();
	for(int i=0;i<100;i++){
		testc[i]->completeOperation();
	}
	log("compress %.4f",get_time_elapsed(sst)/100);

	log("%d vertices %d edges %d faces",compressed->size_of_vertices(), compressed->size_of_halfedges()/2, compressed->true_triangle_size());

	log("start decompressing");

	HiMesh *himesh;
	char path[256];
	int itertime = 100;
	for(int i=start_lod;i<=end_lod;i++){
		int lod = 10*i;
		MyMesh *tested[itertime];
		for(int t=0;t<itertime;t++){
			tested[t] = hispeed::decompress_mesh(compressed, lod, false);
		}
		starttime = get_cur_time();
		for(int t=0;t<itertime;t++){
			tested[t]->completeOperation();
		}
		double testedtime = get_time_elapsed(starttime,true);

		for(int t=0;t<itertime;t++){
			delete tested[t];
		}


		MyMesh *decompressed = hispeed::decompress_mesh(compressed, lod);
		decompressed->completeOperation();
		logt("decompress %3d lod %5d vertices %5d edges %5d faces avg(%.4f)", starttime, lod,
				decompressed->size_of_vertices(), decompressed->size_of_halfedges()/2, decompressed->true_triangle_size(), testedtime/itertime);
		sprintf(path,"/gisdata/lod.%d.off", lod);
		decompressed->writeMeshOff(path);
		if(lod==100){
			float *vertices;
			himesh = new HiMesh(decompressed->p_data,decompressed->dataOffset);
		}
		delete decompressed;
	}
	delete compressed;
	starttime = get_cur_time();
	float *vertices = NULL;
	himesh->advance_to(100);
	logt("decompress", starttime);
	size_t size = himesh->fill_vertices(vertices);
	logt("fill vertices %d with %d bytes (%ld bytes)", starttime, size, size*3*sizeof(float),himesh->dataOffset);
	char *zcomp = new char[size*3*sizeof(float)];
	unsigned long compressedsize;
	for(int i=1;i<10;i++){
		int nResult = compress2((unsigned char *)zcomp, &compressedsize, (unsigned char *)vertices, size*3*sizeof(float),i);
		logt("compress %d level %ld bytes",starttime,i,compressedsize);
	}

	delete []vertices;

}
