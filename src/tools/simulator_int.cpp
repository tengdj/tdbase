/*
 * data_generator.cpp
 *
 *  Created on: Nov 27, 2019
 *      Author: teng
 *
 *  generate testing data by duplicating some given examples
 *
 */

#include <fstream>
#include <tuple>
#include "himesh.h"
#include "tile.h"
#include "popl.h"

using namespace tdbase;
using namespace std;
using namespace popl;

// 50M for each vessel and the nucleis around it
const int buffer_size = 50*(1<<20);
HiMesh *nuclei = NULL;
map<int, HiMesh *> nucleis;

aab nuclei_box;
int num_nuclei = 10000;
int amplify_ratio = 10;

vector<HiMesh_Wrapper *> generated_nucleis;

pthread_mutex_t mylock;

bool multi_lods = false;

void load_prototype(const char *nuclei_path){

	// load the nuclei
	if(multi_lods){
		char path[256];
		// load the nuclei
		for(int lod=100;lod>=20;lod-=20){
			sprintf(path, "%s_%d.off", nuclei_path, lod);
			log("loading %s",path);
			HiMesh *m = read_mesh(path);
			assert(m);
			nucleis[lod] = m;
			if(lod == 100){
				nuclei = m;
			}
		}
	}else{
		nuclei = read_mesh(nuclei_path);
	}
	assert(nuclei);
	auto mbb = nuclei->get_mbb();
	nuclei_box = nuclei->shift(-mbb.low[0], -mbb.low[1], -mbb.low[2]);
}

HiMesh_Wrapper *organize_data(HiMesh *mesh, float shift[3]){

	HiMesh *local_mesh = mesh->clone_mesh();
	local_mesh->shift(shift[0], shift[1], shift[2]);
	local_mesh->encode();

	HiMesh_Wrapper *wr = new HiMesh_Wrapper(local_mesh);
	return wr;
}

HiMesh_Wrapper *organize_data(map<int, HiMesh *> &meshes, float shift[3]){

	map<int, HiMesh *> local_meshes;

	for(auto m:meshes){
		HiMesh *nmesh = m.second->clone_mesh();
		nmesh->shift(shift[0], shift[1], shift[2]);
		local_meshes[m.first] = nmesh;
	}

	HiMesh_Wrapper *wr = new HiMesh_Wrapper(local_meshes);
	return wr;
}

/*
 * generate the binary data of nucleis around vessel with
 * a given shift base
 *
 * */
void* generate_unit(void* arg) {
	log("thread is started");

	bool stop = false;
	float shift[3];
	while(!stop){
		
		shift[0] = tdbase::get_rand_double() * nuclei_box.high[0] * amplify_ratio;
		shift[1] = tdbase::get_rand_double() * nuclei_box.high[1] * amplify_ratio;
		shift[2] = tdbase::get_rand_double() * nuclei_box.high[2] * amplify_ratio;

		HiMesh_Wrapper *wr;
		if(multi_lods){
			wr = organize_data(nucleis, shift);
		}else{
			wr = organize_data(nuclei, shift);
		}

		pthread_mutex_lock(&mylock);
		if (generated_nucleis.size()<num_nuclei) {
			generated_nucleis.push_back(wr);
		}
		else {
			delete wr;
			stop = true;
		}
		pthread_mutex_unlock(&mylock);
	}
	return NULL;
}

int main(int argc, char **argv){

// argument parsing
	string nuclei_pt;
	string output_path;
	int num_threads;


	OptionParser op("Simulator");
	auto help_option = op.add<Switch>("h", "help", "produce help message");
	auto hausdorff_option = op.add<Switch>("", "hausdorff", "enable Hausdorff distance calculation", &HiMesh::use_hausdorff);
	auto multi_lods_option = op.add<Switch>("m", "multi_lods", "the input are polyhedrons in multiple files", &multi_lods);
	op.add<Value<string>>("n", "nuclei", "path to the nuclei prototype file", "nuclei.pt", &nuclei_pt);
	op.add<Value<string>>("o", "output", "prefix of the output files", "default", &output_path);
	op.add<Value<int>>("t", "threads", "number of threads", tdbase::get_num_threads(), &num_threads);
	op.add<Value<int>>("", "amplify_ratio", "how big, in terms of nuclei size, in each dimension", 10, &amplify_ratio);
	op.add<Value<int>>("", "n", "number of nuclei", 10000, &num_nuclei);
	op.add<Value<int>>("", "verbose", "verbose level", 0, &config.verbose);
	op.add<Value<uint32_t>>("r", "sample_rate", "sampling rate for Hausdorff distance calculation", 30, &HiMesh::sampling_rate);
	op.parse(argc, argv);

// processing section

	struct timeval start = get_cur_time();

	pthread_mutex_init(&mylock, NULL);
	char nuclei_output[256];
	char config[100];
	sprintf(config,"nu%d_%d_r%d",
			num_nuclei,amplify_ratio, HiMesh::sampling_rate);

	sprintf(nuclei_output,"%s_n_%s.dt",output_path.c_str(),config);

	load_prototype(nuclei_pt.c_str());
	logt("load prototype files", start);

	pthread_t threads[num_threads];
	for(int i=0;i<num_threads;i++){
		pthread_create(&threads[i], NULL, generate_unit, NULL);
	}
	for(int i = 0; i < num_threads; i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}
	logt("%ld nucleis are generated",start, generated_nucleis.size());

	Tile *nuclei_tile = new Tile(generated_nucleis);

	if(multi_lods){
		nuclei_tile->dump_raw(nuclei_output);
	}else{
		nuclei_tile->dump_compressed(nuclei_output);
	}

	// clear
	delete nuclei;
}

