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
	struct timeval starttime = get_cur_time();
	MyMesh *compressed = read_mesh();
	assert(compressed->size_of_border_edges()&&"must be manifold");
	log("%d vertices %d edges %d faces",compressed->size_of_vertices(), compressed->size_of_halfedges()/2, compressed->size_of_facets());
	compressed->completeOperation();
	logt("compress", starttime);
	log("start decompressing");

	HiMesh *himesh;
	for(int i=start_lod;i<=end_lod;i++){
		int lod = 10*i;
		MyMesh *decompressed = hispeed::decompress_mesh(compressed, lod);
//		sprintf(path,"lod%d.off", lod);
//		decompressed->writeMeshOff(path);
		logt("decompress %3d lod %5d vertices %5d edges %5d faces", starttime, lod,
				decompressed->size_of_vertices(), decompressed->size_of_halfedges()/2, decompressed->size_of_facets());
		if(lod==100){
			float *vertices;
			himesh = new HiMesh(decompressed->p_data,decompressed->dataOffset);
		}
		delete decompressed;
	}
	delete compressed;

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



