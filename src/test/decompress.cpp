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

	for(int i=start_lod;i<=end_lod;i++){
		int lod = 10*i;
		MyMesh *decompressed = hispeed::decompress_mesh(compressed, lod);
//		sprintf(path,"lod%d.off", lod);
//		decompressed->writeMeshOff(path);
		logt("decompress %3d lod %5d vertices %5d edges %5d faces", starttime, lod,
				decompressed->size_of_vertices(), decompressed->size_of_halfedges()/2, decompressed->size_of_facets());
		delete decompressed;
	}

	delete compressed;

}



