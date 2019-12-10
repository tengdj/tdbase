/*
 * decomp.cpp
 *
 *  Created on: Oct 22, 2019
 *      Author: teng
 */


#include "../PPMC/ppmc.h"
#include "../util/util.h"
#include "../spatial/spatial.h"

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
	compressed->completeOperation();
	logt("compress", starttime);
	log("start decompressing");

	char path[256];
	for(int i=start_lod;i<=end_lod;i++){
		int lod = 10*i;
		MyMesh *decompressed = hispeed::decompress_mesh(compressed, lod);
		sprintf(path,"lod%d.off", lod);
		decompressed->writeMeshOff(path);
		delete decompressed;
		logt("decompress %d", starttime, lod);
	}

	delete compressed;

}



