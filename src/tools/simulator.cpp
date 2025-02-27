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
HiMesh *vessel = NULL;
HiMesh *nuclei = NULL;

map<int, HiMesh *> vessels;
map<int, HiMesh *> nucleis;

aab nuclei_box;
aab vessel_box;

bool *vessel_taken;
int total_slots = 0;

int num_nuclei_per_vessel = 100;
int num_vessel = 50;
int voxel_size = 100;

vector<HiMesh_Wrapper *> generated_nucleis;
vector<HiMesh_Wrapper *> generated_vessels;

pthread_mutex_t mylock;

bool multi_lods = false;
bool allow_intersection = false;

void load_prototype(const char *nuclei_path, const char *vessel_path){

	// load the vessel
	if(multi_lods){
		char path[256];
		for(int lod=100;lod>=20;lod-=20){
			sprintf(path, "%s_%d.off", vessel_path, lod);
			log("loading %s",path);
			HiMesh *m = read_mesh(path);
			assert(m);
			vessels[lod] = m;
			if(lod == 100){
				vessel = m;
			}
		}
	}else {
		vessel = read_mesh(vessel_path);
	}
	assert(vessel);
	aab mbb = vessel->get_mbb();
	vessel_box = vessel->shift(-mbb.low[0], -mbb.low[1], -mbb.low[2]);

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
	if (allow_intersection) {
		nuclei->expand(3);
	}
	mbb = nuclei->get_mbb();
	nuclei_box = nuclei->shift(-mbb.low[0], -mbb.low[1], -mbb.low[2]);

	// how many slots in each dimension can one vessel holds
	int nuclei_num[3];
	for(int i=0;i<3;i++){
		nuclei_num[i] = (int)(vessel_box.high[i]/nuclei_box.high[i]);
	}

	total_slots = nuclei_num[0]*nuclei_num[1]*nuclei_num[2];
	vessel_taken = new bool[total_slots];
	for(int i=0;i<total_slots;i++){
		vessel_taken[i] = false;
	}

	// just for mark the voxels that are taken
	vector<Voxel *> vessel_voxels = vessel->generate_voxels_skeleton(50);

	for(Voxel *v:vessel_voxels){
		int xstart = v->low[0]/vessel_box.high[0]*nuclei_num[0];
		int xend = v->high[0]/vessel_box.high[0]*nuclei_num[0];
		int ystart = v->low[1]/vessel_box.high[1]*nuclei_num[1];
		int yend = v->high[1]/vessel_box.high[1]*nuclei_num[1];
		int zstart = v->low[2]/vessel_box.high[2]*nuclei_num[2];
		int zend = v->high[2]/vessel_box.high[2]*nuclei_num[2];

		xend = min(xend, nuclei_num[0]-1);
		yend = min(yend, nuclei_num[1]-1);
		zend = min(zend, nuclei_num[2]-1);

		//log("%d %d, %d %d, %d %d", xstart, xend, ystart, yend, zstart, zend);
		for(int z=zstart;z<=zend;z++){
			for(int y=ystart;y<=yend;y++){
				for(int x=xstart;x<=xend;x++){
					assert(z*nuclei_num[0]*nuclei_num[1]+y*nuclei_num[0]+x < total_slots);
					vessel_taken[z*nuclei_num[0]*nuclei_num[1]+y*nuclei_num[0]+x] = true;
				}
			}
		}
	}
}

HiMesh_Wrapper *organize_data(HiMesh *mesh, float shift[3]){

	HiMesh *local_mesh = mesh->clone_mesh();
	local_mesh->shift(shift[0], shift[1], shift[2]);

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
int generate_nuclei(float base[3]){
	int nuclei_num[3];
	for(int i=0;i<3;i++){
		nuclei_num[i] = (int)(vessel_box.high[i]/nuclei_box.high[i]);
	}
	float shift[3];
	int generated = 0;
	bool *taken = new bool[total_slots];
	memcpy(taken, vessel_taken, total_slots*sizeof(bool));

	int available_count = 0;
	for(int i=0;i<total_slots;i++){
		if(!vessel_taken[i]){
			available_count++;
		}
	}
	//log("%d %d %d",total_slots,num_nuclei_per_vessel,available_count);
	assert(available_count>num_nuclei_per_vessel && "should have enough slots");

	if (!allow_intersection) {
		const int gap = available_count / num_nuclei_per_vessel;

		int skipped = 0;
		int idx = 0;

		while (idx < total_slots) {
			if (taken[idx]) {
				idx++;
				continue;
			}
			if (skipped++ == gap) {
				generated++;
				skipped = 0;
				taken[idx] = true;

				int z = idx / (nuclei_num[0] * nuclei_num[1]);
				int y = (idx % (nuclei_num[0] * nuclei_num[1])) / nuclei_num[0];
				int x = (idx % (nuclei_num[0] * nuclei_num[1])) % nuclei_num[0];

				shift[0] = x * nuclei_box.high[0] + base[0];
				shift[1] = y * nuclei_box.high[1] + base[1];
				shift[2] = z * nuclei_box.high[2] + base[2];

				HiMesh_Wrapper* wr;
				if (multi_lods) {
					wr = organize_data(nucleis, shift);
				}
				else {
					wr = organize_data(nuclei, shift);
				}

				pthread_mutex_lock(&mylock);
				generated_nucleis.push_back(wr);
				pthread_mutex_unlock(&mylock);
			}
			idx++;
		}
	}
	else {
		while (generated++<num_nuclei_per_vessel) {
			shift[0] = tdbase::get_rand_double() * vessel_box.high[0] + base[0];
			shift[1] = tdbase::get_rand_double() * vessel_box.high[1] + base[1];
			shift[2] = tdbase::get_rand_double() * vessel_box.high[2] + base[2];
			HiMesh_Wrapper* wr;
			if (multi_lods) {
				wr = organize_data(nucleis, shift);
			}
			else {
				wr = organize_data(nuclei, shift);
			}
			pthread_mutex_lock(&mylock);
			generated_nucleis.push_back(wr);
			pthread_mutex_unlock(&mylock);
		}
	}
	
	return generated;
}

queue<tuple<float, float, float>> jobs;

long global_generated = 0;

void *generate_unit(void *arg){
	log("thread is started");
	bool complete = false;
	while(true){
		tuple<float, float, float> job;
		pthread_mutex_lock(&mylock);
		complete = jobs.empty();
		if(!complete){
			job = jobs.front();
			jobs.pop();
			log("%ld jobs left", jobs.size());
		}
		pthread_mutex_unlock(&mylock);
		if(complete){
			break;
		}
		size_t offset = 0;
		size_t vessel_offset = 0;
		float base[3] = {get<0>(job), get<1>(job), get<2>(job)};
		int generated = generate_nuclei(base);

		HiMesh_Wrapper *wr;
		if(multi_lods){
			wr = organize_data(vessels, base);
		}else{
			wr = organize_data(vessel, base);
		}
		pthread_mutex_lock(&mylock);
		global_generated += generated;
		generated_vessels.push_back(wr);
		pthread_mutex_unlock(&mylock);
	}
	return NULL;
}

int main(int argc, char **argv){
	string nuclei_pt;
	string vessel_pt;
	string output_path;
	int num_threads;

	pthread_mutex_init(&mylock, NULL);

	OptionParser op("Simulator");
	auto help_option = op.add<Switch>("h", "help", "produce help message");
	auto hausdorff_option = op.add<Switch>("", "hausdorff", "enable Hausdorff distance calculation", &HiMesh::use_hausdorff);
	auto multi_lods_option = op.add<Switch>("m", "multi_lods", "the input are polyhedrons in multiple files", &multi_lods);
	auto allow_intersection_option = op.add<Switch>("i", "allow_intersection", "allow the nuclei to intersect with other nuclei or vessel",&allow_intersection);
	op.add<Value<string>>("n", "nuclei", "path to the nuclei prototype file", "nuclei.pt", &nuclei_pt);
	op.add<Value<string>>("v", "vessel", "path to the vessel prototype file", "vessel.pt", &vessel_pt);
	op.add<Value<string>>("o", "output", "prefix of the output files", "default", &output_path);
	op.add<Value<int>>("t", "threads", "number of threads", tdbase::get_num_threads(), &num_threads);
	op.add<Value<int>>("", "nv", "number of vessels", 50, &num_vessel);
	op.add<Value<int>>("", "nu", "number of nucleis per vessel", 100, &num_nuclei_per_vessel);
	op.add<Value<int>>("", "vs", "number of vertices in each voxel",100, &voxel_size);
	op.add<Value<int>>("", "verbose", "verbose level",0, &global_ctx.verbose);
	op.add<Value<uint32_t>>("r", "sample_rate", "sampling rate for Hausdorff distance calculation", 30, &HiMesh::sampling_rate);

	op.parse(argc, argv);

	struct timeval start = get_cur_time();

	char vessel_output[256];
	char nuclei_output[256];
	char config[100];
	sprintf(config,"nv%d_nu%d_vs%d_r%d",
			num_vessel, num_nuclei_per_vessel, voxel_size,
			HiMesh::sampling_rate);

	sprintf(vessel_output,"%s_v_%s.dt",output_path.c_str(),config);
	sprintf(nuclei_output,"%s_n_%s.dt",output_path.c_str(),config);

	// generate some job for worker to process
	int dim1 = (int)(pow((float)num_vessel, (float)1.0/3)+0.5);
	int dim2 = dim1;
	int dim3 = num_vessel/(dim1*dim2);
	int x_dim,y_dim,z_dim;
	if(vessel_box.high[0]-vessel_box.low[0] > vessel_box.high[1]-vessel_box.low[1] && vessel_box.high[0]-vessel_box.low[0] > vessel_box.high[2]-vessel_box.low[2] ){
		x_dim = min(dim1, dim3);
		y_dim = max(dim1, dim3);
		z_dim = max(dim1, dim3);
	}else if(vessel_box.high[1]-vessel_box.low[1] > vessel_box.high[0]-vessel_box.low[0] && vessel_box.high[1]-vessel_box.low[1] > vessel_box.high[2]-vessel_box.low[2] ){
		y_dim = min(dim1, dim3);
		x_dim = max(dim1, dim3);
		z_dim = max(dim1, dim3);
	}else{
		z_dim = min(dim1, dim3);
		x_dim = max(dim1, dim3);
		y_dim = max(dim1, dim3);
	}
	load_prototype(nuclei_pt.c_str(), vessel_pt.c_str());
	logt("load prototype files", start);

	for(int i=0;i<x_dim;i++){
		for(int j=0;j<y_dim;j++){
			for(int k=0;k<z_dim;k++){
				tuple<float, float, float> tp(i*vessel_box.high[0], j*vessel_box.high[1], k*vessel_box.high[2]);
				jobs.push(tp);
			}
		}
	}
	size_t vessel_num = jobs.size();

	pthread_t threads[num_threads];
	for(int i=0;i<num_threads;i++){
		pthread_create(&threads[i], NULL, generate_unit, NULL);
	}
	for(int i = 0; i < num_threads; i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}
	logt("%ld vessels %ld nucleis are generated",start,vessel_num,global_generated);

	Tile *nuclei_tile = new Tile(generated_nucleis);
	Tile *vessel_tile = new Tile(generated_vessels);

	if(multi_lods){
		nuclei_tile->dump_raw(nuclei_output);
		vessel_tile->dump_raw(vessel_output);
	}else{
		nuclei_tile->dump_compressed(nuclei_output);
		vessel_tile->dump_compressed(vessel_output);
	}

	// clear
	delete vessel;
	delete nuclei;
	delete []vessel_taken;
}

