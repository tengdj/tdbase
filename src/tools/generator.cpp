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
HiMesh *vessel = NULL;
vector<Voxel *> vessel_voxels;
bool *vessel_taken;
int total_slots = 0;

aab nuclei_box;
aab vessel_box;

int num_nuclei_per_vessel = 10000;
int num_vessel = 1;
float shrink = 10;
int voxel_size = 100;
int shifted_range = 30;

void load_prototype(const char *nuclei_path, const char *vessel_path){

	// load the vessel
	vector<HiMesh *> vessels = read_meshes(vessel_path, 1);
	assert(vessels.size()==1);
	vessel = vessels[0];
	aab mbb = vessel->get_mbb();
	vessel_box = vessel->shift(-mbb.low[0], -mbb.low[1], -mbb.low[2]);

	// load the nucleis
	nucleis = read_meshes(nuclei_path, 1);
	for(HiMesh *mesh:nucleis){
		mbb = mesh->shrink(shrink);
		mbb = mesh->shift(-mbb.low[0], -mbb.low[1], -mbb.low[2]);
		nuclei_box.update(mbb);
		vector<Voxel *> vxls = mesh->generate_voxels_skeleton(mesh->size_of_vertices()/voxel_size);
		nucleis_voxels.push_back(vxls);
	}

	vessel_box.print();
	nuclei_box.print();

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
	vector<Voxel *> vessel_voxels = vessel->generate_voxels_skeleton(1000);;
	hispeed::write_voxels(vessel_voxels, "/gisdata/3dpro/generated/voxels.OFF");
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
void organize_data(HiMesh *mesh, vector<Voxel *> &voxels, float shift[3], char *data, size_t &offset){
	HiMesh *nmesh = mesh->clone_mesh();
	nmesh->shift(shift[0], shift[1], shift[2]);
	nmesh->encode();
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
	assert(offset<buffer_size);
	delete nmesh;
}

/*
 * generate the binary data of nucleis around vessel with
 * a given shift base
 *
 * */
int generate_nuclei(float base[3], char *data, size_t &offset, char *data2, size_t &offset2){
	int nuclei_num[3];
	for(int i=0;i<3;i++){
		nuclei_num[i] = (int)(vessel_box.high[i]/nuclei_box.high[i]);
	}
	float shift[3];
	int generated = 0;
	bool *taken = new bool[total_slots];
	memcpy(taken, vessel_taken, total_slots*sizeof(bool));

	int taken_count = 0;
	for(int i=0;i<total_slots;i++){
		if(vessel_taken[i]){
			taken_count++;
		}
	}
	//log("%d %d %d",total_slots,num_nuclei_per_vessel,taken_count);
	assert(total_slots>num_nuclei_per_vessel+taken_count && "should have enough slots");
	while(generated++<num_nuclei_per_vessel){
		int idx = get_rand_number(total_slots);
		int oidx = idx;
		for(;taken[idx];idx = (idx+1)%total_slots);
		taken[idx] = true;

		int z = idx/(nuclei_num[0]*nuclei_num[1]);
		int y = (idx%(nuclei_num[0]*nuclei_num[1]))/nuclei_num[0];
		int x = (idx%(nuclei_num[0]*nuclei_num[1]))%nuclei_num[0];

		shift[0] = x*nuclei_box.high[0]+base[0];
		shift[1] = y*nuclei_box.high[1]+base[1];
		shift[2] = z*nuclei_box.high[2]+base[2];

		int polyid = get_rand_number(nucleis.size()-1);
		organize_data(nucleis[polyid], nucleis_voxels[polyid], shift, data, offset);

		float shift2[3];
		shift2[0] = shift[0]+nuclei_box.high[0]*(get_rand_number(shifted_range)*1.0)/100.0*(flip_coin()?1:-1);
		shift2[1] = shift[1]+nuclei_box.high[1]*(get_rand_number(shifted_range)*1.0)/100.0*(flip_coin()?1:-1);
		shift2[2] = shift[2]+nuclei_box.high[2]*(get_rand_number(shifted_range)*1.0)/100.0*(flip_coin()?1:-1);

		int polyid2 = get_rand_number(nucleis.size()-1);
		organize_data(nucleis[polyid2], nucleis_voxels[polyid2], shift2, data2, offset2);
		//log("%d %d %d %d",idx,oidx,polyid,polyid2);
	}
	return generated;
}

ofstream *v_os = NULL;
ofstream *os = NULL;
ofstream *os2 = NULL;

queue<tuple<float, float, float>> jobs;
pthread_mutex_t mylock;

long global_generated = 0;

void *generate_unit(void *arg){
	char *vessel_data = new char[buffer_size];
	char *data = new char[buffer_size];
	char *data2 = new char[buffer_size];
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
		size_t offset2 = 0;
		size_t vessel_offset = 0;
		float base[3] = {get<0>(job), get<1>(job), get<2>(job)};
		int generated = generate_nuclei(base, data, offset, data2, offset2);
		organize_data(vessel, vessel_voxels, base, vessel_data, vessel_offset);
		pthread_mutex_lock(&mylock);
		os->write(data, offset);
		os2->write(data2, offset2);
		v_os->write(vessel_data, vessel_offset);
		global_generated += generated;
		pthread_mutex_unlock(&mylock);
	}
	delete vessel_data;
	delete data;
	delete data2;
	return NULL;
}

int main(int argc, char **argv){
	string nuclei_pt = "../data/nuclei.pt";
	string vessel_pt = "../data/vessel.pt";
	string output_path;
	int num_threads = hispeed::get_num_threads();

	pthread_mutex_init(&mylock, NULL);

	po::options_description desc("joiner usage");
	desc.add_options()
		("help,h", "produce help message")
		("nuclei,u", po::value<string>(&nuclei_pt), "path to the nuclei prototype file")
		("vessel,v", po::value<string>(&vessel_pt), "path to the vessel prototype file")
		("output,o", po::value<string>(&output_path)->required(), "prefix of the output files")
		("shrink,s", po::value<float>(&shrink), "shrink the size of nuclei by how many times")
		("threads,n", po::value<int>(&num_threads), "number of threads")
		("shifted_range,r", po::value<int>(&shifted_range), "range of the second nucleus can be shifted")
		("nv", po::value<int>(&num_vessel), "number of vessels")
		("nu", po::value<int>(&num_nuclei_per_vessel), "number of nucleis per vessel")
		("vs", po::value<int>(&voxel_size), "number of vertices in each voxel")
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
	char nuclei_output2[256];
	char config[100];
	sprintf(config,"nv%d_nu%d_s%d_vs%d_sr%d_r%d_cm%d",
			num_vessel, num_nuclei_per_vessel,
			(int)shrink, voxel_size, shifted_range,
			HiMesh::sampling_rate, HiMesh::calculate_method);

	sprintf(vessel_output,"%s_v_%s.mt",output_path.c_str(),config);
	remove(vessel_output);
	sprintf(nuclei_output,"%s_n_%s.mt",output_path.c_str(),config);
	remove(nuclei_output);
	sprintf(nuclei_output2,"%s_n2_%s.mt",output_path.c_str(),config);
	remove(nuclei_output2);

	sprintf(vessel_output,"%s_v_%s.dt",output_path.c_str(),config);
	sprintf(nuclei_output,"%s_n_%s.dt",output_path.c_str(),config);
	sprintf(nuclei_output2,"%s_n2_%s.dt",output_path.c_str(),config);

	v_os = new std::ofstream(vessel_output, std::ios::out | std::ios::binary);
	os = new std::ofstream(nuclei_output, std::ios::out | std::ios::binary);
	os2 = new std::ofstream(nuclei_output2, std::ios::out | std::ios::binary);

	// generate some job for worker to process
	int x_dim = (int)pow((float)num_vessel, 1.0/3);
	int y_dim = x_dim;
	int z_dim = num_vessel/(x_dim*y_dim);
	num_nuclei_per_vessel = (num_nuclei_per_vessel*num_vessel)/(x_dim*y_dim*z_dim);
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

	// close the output stream
	v_os->close();
	os->close();
	os2->close();
	delete v_os;
	delete os;
	delete os2;

	// clear the vessel related objects
	delete vessel;
	for(Voxel *v:vessel_voxels){
		delete v;
	}
	vessel_voxels.clear();

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
	delete []vessel_taken;
}

