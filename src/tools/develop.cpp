/*
 * develop.cpp
 *
 *  Created on: Nov 17, 2022
 *      Author: teng
 */

#include "../mymesh/mymesh.h"

using namespace tmesh;

// parse the metadata
char *load_one(const char *path){
	FILE *dt_fs = fopen(path, "r");
	assert(dt_fs);
	size_t dsize = 0;
	long offset = 0;
	size_t index = 0;
	fseek(dt_fs, 0, SEEK_SET);
	assert(fread((void *)&dsize, sizeof(size_t), 1, dt_fs));
	char *data = new char[dsize];
	assert(fread(data, sizeof(char), dsize, dt_fs) == dsize);
	fclose(dt_fs);
	return data;
}

int main(int argc, char **argv){

	char *dat = load_one("/home/teng/git/3DPro/src/barycenter_n_nv50_nu200_s10_vs100_r30.dt");
	tmesh::TMesh *poly = new tmesh::TMesh(0);
	poly->load(dat, true);

//	poly->load(argv[1]);
//	poly->remove_orphan_vertices();
//
//	char path[256];
//	sprintf(path,"/gisdata/mymesh/original.off");
//	poly->dumpto(path);
//	poly->evaluate();
//
//	for(int i=0;i<10;i++){
//		poly->reset_states();
//		struct timeval start = get_cur_time();
//		poly->compress();
//		logt("compress %d", start, i);
//		//poly->print();
//		sprintf(path,"/gisdata/mymesh/compressed_%d.off", i);
//		poly->dumpto(path);
//		poly->evaluate();
//	}
//	delete poly;

	return 0;
}





