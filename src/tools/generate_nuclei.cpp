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

using namespace hispeed;
using namespace std;

namespace po = boost::program_options;

// 50M for each vessel and the nucleis around it
const int buffer_size = 50*(1<<20);
vector<HiMesh *> nucleis;
vector<vector<Voxel *>> nucleis_voxels;
aab nuclei_box;
aab global_space;

void load_prototype(const char *nuclei_path){

	// load the nucleis
	nucleis = read_meshes(nuclei_path, 100);
	for(HiMesh *mesh:nucleis){
		aab mbb = mesh->get_mbb();
		mbb = mesh->shift(-mbb.low[0], -mbb.low[1], -mbb.low[2]);
		nuclei_box.update(mbb);
		vector<Voxel *> vxls = mesh->generate_voxels_skeleton(mesh->size_of_vertices()/500);
		nucleis_voxels.push_back(vxls);
	}
	nuclei_box.print();
}

/*
 * generate the binary data for a polyhedron and its voxels
 * */
int idx = 0;
void organize_data(HiMesh *mesh, vector<Voxel *> &voxels, float shift[3], char *data, size_t &offset){
	HiMesh *nmesh = mesh->clone_mesh();
	nmesh->shift(shift[0], shift[1], shift[2]);
	nmesh->encode();

	process_lock();
	char path[256];
	sprintf(path, "/gisdata/nuclei_%d.off", idx++);
	nmesh->write_to_off(path);
	//hispeed::write_polyhedron(&shifted, ids++);
	size_t size = nmesh->get_data_size();
	memcpy(data+offset, (char *)&size, sizeof(size_t));
	offset += sizeof(size_t);
	memcpy(data+offset, nmesh->get_data(), nmesh->get_data_size());
	offset += nmesh->get_data_size();
	size = voxels.size();
	memcpy(data+offset, (char *)&size, sizeof(size_t));
	offset += sizeof(size_t);
	float box_tmp[3];
	for(Voxel *v:voxels){
		for(int i=0;i<3;i++){
			box_tmp[i] = v->low[i]+shift[i];
		}
		memcpy(data+offset, (char *)box_tmp, 3*sizeof(float));
		offset += 3*sizeof(float);

		for(int i=0;i<3;i++){
			box_tmp[i] = v->high[i]+shift[i];
		}
		memcpy(data+offset, (char *)box_tmp, 3*sizeof(float));
		offset += 3*sizeof(float);

		for(int i=0;i<3;i++){
			box_tmp[i] = v->core[i]+shift[i];
		}
		memcpy(data+offset, (char *)box_tmp, 3*sizeof(float));
		offset += 3*sizeof(float);
	}
	process_unlock();
	assert(offset<buffer_size);
	delete nmesh;
}

ofstream *os = NULL;

queue<int> jobs;
pthread_mutex_t mylock;

void *generate_unit(void *arg){
	char *data = new char[buffer_size];
	size_t offset = 0;
	bool complete = false;
	int generated = 0;
	while(true){
		int job;
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

		float shift[3];
		shift[0] = hispeed::get_rand_double()*global_space.high[0];
		shift[1] = hispeed::get_rand_double()*global_space.high[1];
		shift[2] = hispeed::get_rand_double()*global_space.high[2];

		int polyid = get_rand_number(nucleis.size()-1);
		organize_data(nucleis[polyid], nucleis_voxels[polyid], shift, data, offset);

		if(++generated > 10){
			pthread_mutex_lock(&mylock);
			os->write(data, offset);
			pthread_mutex_unlock(&mylock);
			offset = 0;
			generated = 0;
		}
	}

	if(offset>0){
		pthread_mutex_lock(&mylock);
		os->write(data, offset);
		pthread_mutex_unlock(&mylock);
	}

	delete data;
	return NULL;
}

int main(int argc, char **argv){
	string nuclei_pt = "../data/nuclei.pt";
	string output_path;
	int num_threads = hispeed::get_num_threads();
	int num_objects = 10000;
	int density = 10;

	pthread_mutex_init(&mylock, NULL);

	po::options_description desc("joiner usage");
	desc.add_options()
		("help,h", "produce help message")
		("nuclei,u", po::value<string>(&nuclei_pt), "path to the nuclei prototype file")
		("output,o", po::value<string>(&output_path)->required(), "prefix of the output files")
		("threads,n", po::value<int>(&num_threads), "number of threads")
		("num_objects", po::value<int>(&num_objects), "number of objects")
		("density,d", po::value<int>(&density), "the density of objects (the gap between objects)")
		("verbose", po::value<int>(&global_ctx.verbose), "verbose level")
		("sample_rate", po::value<uint>(&HiMesh::sampling_rate), "sampling rate for Hausdorff distance calculation")
		("calculate_method", po::value<int>(&HiMesh::calculate_method), "hausdorff distance calculating method [0ALL|1BVH|2ASSOCIATE|3ASSOCIATE_CYLINDER|4NULL]")
		;

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
	sprintf(config,"no%d_ds%d_r%d_cm%d",
			num_objects, density,
			HiMesh::sampling_rate, HiMesh::calculate_method);

	int eachdim = pow(num_objects, 1.0/3)+1;

	num_objects = eachdim*eachdim*eachdim;

	sprintf(nuclei_output,"%s_%s.mt",output_path.c_str(),config);
	remove(nuclei_output);
	sprintf(nuclei_output,"%s_%s.dt",output_path.c_str(),config);

	os = new std::ofstream(nuclei_output, std::ios::out | std::ios::binary);

	load_prototype(nuclei_pt.c_str());
	logt("load prototype files", start);

	for(int i=0;i<3;i++){
		global_space.low[i] = 0;
		global_space.high[i] = nuclei_box.high[i]*eachdim*(density+1);
	}

	for(int i=0;i<num_objects;i++){
		jobs.push(i);
	}

	pthread_t threads[num_threads];
	for(int i=0;i<num_threads;i++){
		pthread_create(&threads[i], NULL, generate_unit, NULL);
	}
	for(int i = 0; i < num_threads; i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}
	logt("%ld nucleis are generated",start,num_objects);

	// close the output stream
	os->close();
	delete os;
	// clear the nuclei related objects
	for(HiMesh *n:nucleis){
		delete n;
	}
	nucleis.clear();
	for(vector<Voxel *> &vs:nucleis_voxels){
		for(Voxel *v:vs){
			delete v;
		}
		vs.clear();
	}
	nucleis_voxels.clear();
}

