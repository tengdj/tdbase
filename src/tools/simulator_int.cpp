/*
 * data_generator.cpp
 *
 *  Created on: Nov 27, 2019
 *      Author: teng
 *
 *  generate testing data by duplicating some given examples
 *
 */

#include <boost/program_options.hpp>
#include <fstream>
#include <tuple>
#include "himesh.h"
#include "tile.h"

using namespace tdbase;
using namespace std;

namespace po = boost::program_options;

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
	string nuclei_pt = "../data/nuclei.pt";
	string output_path;
	int num_threads = tdbase::get_num_threads();

	pthread_mutex_init(&mylock, NULL);

	int cm = HiMesh::calculate_method;
	po::options_description desc("joiner usage");
	desc.add_options()
		("help,h", "produce help message")
		("ppvp,p", "enable the ppvp mode, for simulator and join query")
		("multi_lods,m", "the input are polyhedrons in multiple files")
		("nuclei,u", po::value<string>(&nuclei_pt), "path to the nuclei prototype file")
		("output,o", po::value<string>(&output_path)->required(), "prefix of the output files")
		("threads,n", po::value<int>(&num_threads), "number of threads")
		("nu", po::value<int>(&num_nuclei), "number of nuclei")
		("amplify_ratio", po::value<int>(&amplify_ratio), "how big, in terms of nuclei size, in each dimension")
		("verbose", po::value<int>(&global_ctx.verbose), "verbose level")
		("sample_rate,r", po::value<uint32_t>(&HiMesh::sampling_rate), "sampling rate for Hausdorff distance calculation (default 30)")
		("calculate_method", po::value<int>(&cm), "hausdorff distance calculating method [0NULL|1BVH(default)|2ASSOCIATE|3ASSOCIATE_CYLINDER]")
		;
	HiMesh::calculate_method = (Hausdorff_Computing_Type)cm;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);

	if(vm.count("help")){
		cout<<desc<<endl;
		return 0;
	}
	po::notify(vm);

	if(vm.count("ppvp")){
		global_ctx.ppvp = true;
	}
	if(vm.count("multi_lods")){
		multi_lods = true;
	}

	struct timeval start = get_cur_time();

	char nuclei_output[256];
	char config[100];
	sprintf(config,"nu%d_%d_r%d_cm%d",
			num_nuclei,amplify_ratio, HiMesh::sampling_rate, HiMesh::calculate_method);

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

