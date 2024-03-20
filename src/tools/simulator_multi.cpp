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
#include "tile.h"
#include "himesh.h"

using namespace hispeed;
using namespace std;

namespace po = boost::program_options;

// 50M for each vessel and the nucleis around it
const int buffer_size = 50*(1<<20);
HiMesh *vessel = NULL;
HiMesh *nuclei = NULL;

map<int, HiMesh *> vessels;
map<int, HiMesh *> nucleis;

aab nuclei_box;
aab vessel_box;

vector<Voxel *> vessel_voxels;
vector<Voxel *> nuclei_voxels;

bool *vessel_taken;
int total_slots = 0;

int num_nuclei_per_vessel = 100;
int num_vessel = 50;
int voxel_size = 100;

vector<HiMesh_Wrapper *> generated_nucleis;
vector<HiMesh_Wrapper *> generated_vessels;

void load_prototype(const char *nuclei_path, const char *vessel_path){

	// load the vessel
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

	aab mbb = vessel->get_mbb();
	for(auto v:vessels){
		v.second->shift(-mbb.low[0], -mbb.low[1], -mbb.low[2]);
	}
	vessel_box = vessel->get_mbb();

	mbb = nuclei->get_mbb();
	for(auto n:nucleis){
		n.second->shift(-mbb.low[0], -mbb.low[1], -mbb.low[2]);
	}
	nuclei_box = nuclei->get_mbb();

	nuclei_voxels = nuclei->generate_voxels_skeleton(nuclei->size_of_vertices()/voxel_size);

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
	vector<Voxel *> vessel_voxels = vessel->generate_voxels_skeleton(1000);

	//hispeed::write_voxels(vessel_voxels, "/gisdata/3dpro/generated/voxels.OFF");
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

/*
 * generate the binary data for a polyhedron and its voxels
 * */
HiMesh_Wrapper *organize_data(map<int, HiMesh *> &meshes, vector<Voxel *> &voxels, float shift[3]){

	map<int, HiMesh *> local_meshes;
	vector<Voxel *> local_voxels;

	for(auto m:meshes){
		HiMesh *nmesh = m.second->clone_mesh();
		nmesh->shift(shift[0], shift[1], shift[2]);
		local_meshes[m.first] = nmesh;
	}

	for(Voxel *v:voxels){
		Voxel *vv = new Voxel();
		for(int i=0;i<3;i++){
			vv->low[i] = v->low[i]+shift[i];
		}
		for(int i=0;i<3;i++){
			vv->high[i] = v->high[i]+shift[i];
		}
		for(int i=0;i<3;i++){
			vv->core[i] = v->core[i]+shift[i];
		}
		local_voxels.push_back(vv);
	}

	HiMesh_Wrapper *wr = new HiMesh_Wrapper(local_meshes, local_voxels);
	return wr;
}
pthread_mutex_t mylock;

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
	//log("%d %d %d",total_slots,num_nuclei_per_vessel,taken_count);
	assert(available_count<num_nuclei_per_vessel && "should have enough slots");


	const int gap = available_count/num_nuclei_per_vessel;

	int skipped = 0;
	int idx = 0;

	while(idx<total_slots){
		if(taken[idx]){
			idx++;
			continue;
		}
		if(skipped++ == gap){
			generated++;
			skipped = 0;
			taken[idx] = true;

			int z = idx/(nuclei_num[0]*nuclei_num[1]);
			int y = (idx%(nuclei_num[0]*nuclei_num[1]))/nuclei_num[0];
			int x = (idx%(nuclei_num[0]*nuclei_num[1]))%nuclei_num[0];

			shift[0] = x*nuclei_box.high[0]+base[0];
			shift[1] = y*nuclei_box.high[1]+base[1];
			shift[2] = z*nuclei_box.high[2]+base[2];

			HiMesh_Wrapper *wr = organize_data(nucleis, nuclei_voxels, shift);
			pthread_mutex_lock(&mylock);
			generated_nucleis.push_back(wr);
			pthread_mutex_unlock(&mylock);
		}
		idx++;
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
		float base[3] = {get<0>(job), get<1>(job), get<2>(job)};
		int generated = generate_nuclei(base);
		HiMesh_Wrapper *hw = organize_data(vessels, vessel_voxels, base);
		pthread_mutex_lock(&mylock);
		generated_vessels.push_back(hw);
		global_generated += generated;
		pthread_mutex_unlock(&mylock);
	}
	return NULL;
}

int main(int argc, char **argv){
	string nuclei_pt = "../data/nuclei.pt";
	string vessel_pt = "../data/vessel.pt";
	string output_path;
	int num_threads = hispeed::get_num_threads();

	pthread_mutex_init(&mylock, NULL);

	int calcm = HiMesh::calculate_method;
	po::options_description desc("joiner usage");
	desc.add_options()
		("help,h", "produce help message")
		("nuclei,u", po::value<string>(&nuclei_pt), "path to the nuclei prototype file")
		("vessel,v", po::value<string>(&vessel_pt), "path to the vessel prototype file")
		("output,o", po::value<string>(&output_path)->required(), "prefix of the output files")
		("threads,n", po::value<int>(&num_threads), "number of threads")
		("nv", po::value<int>(&num_vessel), "number of vessels")
		("nu", po::value<int>(&num_nuclei_per_vessel), "number of nucleis per vessel")
		("vs", po::value<int>(&voxel_size), "number of vertices in each voxel")
		("verbose", po::value<int>(&global_ctx.verbose), "verbose level")
		("sample_rate,r", po::value<uint32_t>(&HiMesh::sampling_rate), "sampling rate for Hausdorff distance calculation (default 30)")
		("calculate_method,m", po::value<int>(&calcm), "hausdorff distance calculating method [0NULL|1BVH(default)|2ASSOCIATE|3ASSOCIATE_CYLINDER]")
		;
	HiMesh::calculate_method = (Hausdorff_Computing_Type)calcm;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);

	if(vm.count("help")){
		cout<<desc<<endl;
		return 0;
	}
	po::notify(vm);

	struct timeval start = get_cur_time();

	char vessel_output[256];
	char nuclei_output[256];
	char config[100];
	sprintf(config,"nv%d_nu%d_vs%d_r%d_cm%d",
			num_vessel, num_nuclei_per_vessel, voxel_size,
			HiMesh::sampling_rate, HiMesh::calculate_method);

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
	vessel_voxels = vessel->generate_voxels_skeleton(voxel_size);

	pthread_t threads[num_threads];
	for(int i=0;i<num_threads;i++){
		pthread_create(&threads[i], NULL, generate_unit, NULL);
	}
	for(int i = 0; i < num_threads; i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}
	logt("%ld vessels with %d voxels %ld nucleis are generated",start,vessel_num,vessel_voxels.size(),global_generated);

	// dump the generated meshes to
	Tile *nuclei_tile = new Tile(generated_nucleis);
	Tile *vessel_tile = new Tile(generated_vessels);

	nuclei_tile->dump_raw(nuclei_output);
	vessel_tile->dump_raw(vessel_output);

	// clear the vessel related objects
	delete vessel;
	for(Voxel *v:vessel_voxels){
		delete v;
	}
	vessel_voxels.clear();

	// clear the nuclei related objects
	delete nuclei;
	for(Voxel *v:nuclei_voxels){
		delete v;
	}
	nuclei_voxels.clear();

	delete []vessel_taken;
}

