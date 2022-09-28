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
vector<Polyhedron *> nucleis;
vector<vector<Voxel *>> nucleis_voxels;
bool *vessel_taken;
Polyhedron *vessel = NULL;
aab nuclei_box;
aab vessel_box;


int num_nuclei_per_vessel = 10000;
int num_vessel = 1;
float shrink = 20;
int voxel_size = 400000;
int shifted_range = 100;

HiMesh *poly_to_himesh(Polyhedron &poly){
	stringstream ss;
	ss<<poly;
	MyMesh *mesh = hispeed::parse_mesh(ss.str(), true);
	HiMesh *himesh = new HiMesh(mesh->p_data, mesh->dataOffset);
	himesh->advance_to(100);
	delete mesh;
	return himesh;
}

MyMesh *poly_to_mesh(Polyhedron &poly){
	stringstream ss;
	ss<<poly;
	return parse_mesh(ss.str(), true);
}

void load_prototype(const char *nuclei_path, const char *vessel_path){

	vector<Voxel *> vessel_voxels;
	// load the vessel
	vector<Polyhedron *> vessels = read_polyhedrons(vessel_path, 1);
	assert(vessels.size()==1);
	vessel = vessels[0];
	vessels.clear();
	aab tmpb;
	for(Polyhedron::Vertex_iterator vi=vessel->vertices_begin();vi!=vessel->vertices_end();vi++){
		Point p = vi->point();
		tmpb.update(p[0], p[1], p[2]);
	}
	for(Polyhedron::Vertex_iterator vi=vessel->vertices_begin();vi!=vessel->vertices_end();vi++){
		Point p = vi->point();
		vi->point() = Point(p[0]-tmpb.low[0], p[1]-tmpb.low[1], p[2]-tmpb.low[2]);
	}
	tmpb.high[0] -= tmpb.low[0];
	tmpb.high[1] -= tmpb.low[1];
	tmpb.high[2] -= tmpb.low[2];
	tmpb.low[0] = 0;
	tmpb.low[1] = 0;
	tmpb.low[2] = 0;
	vessel_box.update(tmpb);
	HiMesh *himesh = poly_to_himesh(*vessel);
	// just for assigning nuclei in the sub space around the vessel
	vessel_voxels = himesh->generate_voxels_skeleton(1);
	delete himesh;

	// load the nucleis
	nucleis = read_polyhedrons(nuclei_path);
	for(Polyhedron *poly:nucleis){
		aab tmpb;
		for(Polyhedron::Vertex_iterator vi=poly->vertices_begin();vi!=poly->vertices_end();vi++){
			Point p = vi->point();
			Point np(p[0]/shrink, p[1]/shrink, p[2]/shrink);
			vi->point() = np;
			tmpb.update(np[0], np[1], np[2]);
		}
		for(Polyhedron::Vertex_iterator vi=poly->vertices_begin();vi!=poly->vertices_end();vi++){
			Point p = vi->point();
			vi->point() = Point((p[0]-tmpb.low[0]), (p[1]-tmpb.low[1]), p[2]-tmpb.low[2]);
		}
		tmpb.high[0] -= tmpb.low[0];
		tmpb.high[1] -= tmpb.low[1];
		tmpb.high[2] -= tmpb.low[2];
		tmpb.low[2] = 0;
		tmpb.low[1] = 0;
		tmpb.low[2] = 0;
		nuclei_box.update(tmpb);
		HiMesh *himesh = poly_to_himesh(*poly);
		vector<Voxel *> vxls = himesh->generate_voxels_skeleton(himesh->size_of_vertices()/voxel_size);
		nucleis_voxels.push_back(vxls);
		delete himesh;
	}

	int nuclei_num[3];
	for(int i=0;i<3;i++){
		nuclei_num[i] = (int)(vessel_box.high[i]/nuclei_box.high[i]);
	}

	int total_slots = nuclei_num[0]*nuclei_num[1]*nuclei_num[2];
	vessel_taken = new bool[total_slots];
	for(Voxel *v:vessel_voxels){
		int xstart = v->low[0]/nuclei_num[0];
		int xend = v->high[0]/nuclei_num[0];
		int ystart = v->low[1]/nuclei_num[1];
		int yend = v->high[1]/nuclei_num[1];
		int zstart = v->low[2]/nuclei_num[2];
		int zend = v->high[2]/nuclei_num[2];
		for(int z=zstart;z<=zend;z++){
			for(int y=ystart;y<=yend;y++){
				for(int x=xstart;x<=xend;x++){
					vessel_taken[z*nuclei_num[0]*nuclei_num[1]+y*nuclei_num[0]+x] = true;
				}
			}
		}
	}
}

Polyhedron shift_polyhedron(float shift[3], Polyhedron &poly_o){
	Polyhedron poly;
	stringstream ss;
	ss << poly_o;
	ss >> poly;
	for(Polyhedron::Vertex_iterator vi=poly.vertices_begin();vi!=poly.vertices_end();vi++){
		Point p = vi->point();
		vi->point() = Point((p[0]+shift[0]), (p[1]+shift[1]), p[2]+shift[2]);
	}
	return poly;
}

/*
 * generate the binary data for a polyhedron and its voxels
 * */

int ids = 0;
inline void organize_data(Polyhedron &poly, vector<Voxel *> voxels,
		float shift[3], char *data, size_t &offset){
	Polyhedron shifted = shift_polyhedron(shift, poly);
	//hispeed::write_polyhedron(&shifted, ids++);
	MyMesh *mesh = poly_to_mesh(shifted);
	memcpy(data+offset, (char *)&mesh->dataOffset, sizeof(size_t));
	offset += sizeof(size_t);
	memcpy(data+offset, mesh->p_data, mesh->dataOffset);
	offset += mesh->dataOffset;
	size_t size = voxels.size();
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
	delete mesh;
}

/*
 * generate the binary data of nucleis around vessel with
 * a given shift base
 *
 * */
inline int generate_nuclei(float base[3], char *data, size_t &offset, char *data2, size_t &offset2){
	int nuclei_num[3];
	for(int i=0;i<3;i++){
		nuclei_num[i] = (int)(vessel_box.high[i]/nuclei_box.high[i]);
	}
	float shift[3];
	int generated = 0;
	int total_slots = nuclei_num[0]*nuclei_num[1]*nuclei_num[2];
	bool *taken = new bool[total_slots];

	assert(total_slots>num_nuclei_per_vessel);

	while(++generated<num_nuclei_per_vessel){
		int idx = get_rand_number(total_slots);
		while(taken[idx]||vessel_taken[idx]){
			idx = get_rand_number(total_slots);
		}
		taken[idx] = true;

		int z = idx/(nuclei_num[0]*nuclei_num[1]);
		int y = (idx%(nuclei_num[0]*nuclei_num[1]))/nuclei_num[0];
		int x = (idx%(nuclei_num[0]*nuclei_num[1]))%nuclei_num[0];

		shift[0] = x*nuclei_box.high[0]+base[0];
		shift[1] = y*nuclei_box.high[1]+base[1];
		shift[2] = z*nuclei_box.high[2]+base[2];

		int polyid = hispeed::get_rand_number(nucleis.size()-1);
		organize_data(*nucleis[polyid], nucleis_voxels[polyid], shift, data, offset);
		{
			float shift2[3];
			shift2[0] = shift[0]+nuclei_box.high[0]*(hispeed::get_rand_number(shifted_range)*1.0)/100.0*(hispeed::get_rand_sample(50)?1:-1);
			shift2[1] = shift[1]+nuclei_box.high[1]*(hispeed::get_rand_number(shifted_range)*1.0)/100.0*(hispeed::get_rand_sample(50)?1:-1);
			shift2[2] = shift[2]+nuclei_box.high[2]*(hispeed::get_rand_number(shifted_range)*1.0)/100.0*(hispeed::get_rand_sample(50)?1:-1);

			int polyid2 = hispeed::get_rand_number(nucleis.size()-1);
			organize_data(*nucleis[polyid2], nucleis_voxels[polyid2], shift2, data2, offset2);
		}
	}
	return generated;
}

ofstream *os = NULL;
ofstream *os2 = NULL;
queue<tuple<float, float, float>> jobs;
pthread_mutex_t mylock;
long global_generated = 0;


void *generate_unit(void *arg){
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
		float base[3] = {get<0>(job),get<1>(job),get<2>(job)};
		int generated = generate_nuclei(base, data, offset, data2, offset2);

		pthread_mutex_lock(&mylock);
		os->write(data, offset);
		os2->write(data2, offset2);
		global_generated += generated;
		pthread_mutex_unlock(&mylock);
	}
	delete data;
	delete data2;
	return NULL;
}

void generate_vessel(const char *path, vector<tuple<float, float, float>> &vessel_shifts){
	struct timeval start = get_cur_time();
	char *data = new char[vessel_shifts.size()*100000*2];
	size_t offset = 0;
	HiMesh *himesh = poly_to_himesh(*vessel);
	vector<Voxel *> voxels = himesh->generate_voxels_skeleton(voxel_size);
	for(tuple<float, float, float> tp:vessel_shifts){
		float shift[3] = {get<0>(tp),get<1>(tp),get<2>(tp)};
		organize_data(*vessel, voxels, shift, data, offset);
	}
	ofstream *v_os = new std::ofstream(path, std::ios::out | std::ios::binary);
	v_os->write(data, offset);
	v_os->close();
	delete v_os;
	logt("%ld vessels are generated with %d voxels",start,vessel_shifts.size(),voxels.size());
	for(Voxel *v:voxels){
		delete v;
	}
	voxels.clear();
	delete himesh;
}


int main(int argc, char **argv){
	string nuclei_pt;
	string vessel_pt;
	string output_path;
	int num_threads = hispeed::get_num_threads();

	po::options_description desc("joiner usage");
	desc.add_options()
		("help,h", "produce help message")
		("nuclei,u", po::value<string>(&nuclei_pt)->required(), "path to the nuclei prototype file")
		("vessel,v", po::value<string>(&vessel_pt)->required(), "path to the vessel prototype file")
		("output,o", po::value<string>(&output_path)->required(), "initial of the output files")
		("shrink,s", po::value<float>(&shrink), "shrink the size of nuclei by how many times")
		("threads,n", po::value<int>(&num_threads), "number of threads")
		("shifted_range,r", po::value<int>(&shifted_range), "range of the second nucleus can be shifted")
		("nv", po::value<int>(&num_vessel), "number of vessels")
		("nu", po::value<int>(&num_nuclei_per_vessel), "number of nucleis per vessel")
		("vs", po::value<int>(&voxel_size), "number of vertices in each voxel")
		;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);

	if(vm.count("help")){
		cout<<desc<<endl;
		return 0;
	}
	po::notify(vm);

	struct timeval start = get_cur_time();

	char nuclei_output[256];
	char nuclei_output2[256];

	char vessel_output[256];
	sprintf(nuclei_output,"%s_n_nv%d_nu%d_s%d_vs%d_r%d.dt",output_path.c_str(),num_vessel,num_nuclei_per_vessel,(int)shrink, voxel_size,shifted_range);
	sprintf(nuclei_output2,"%s_n2_nv%d_nu%d_s%d_vs%d_r%d.dt",output_path.c_str(),num_vessel,num_nuclei_per_vessel,(int)shrink, voxel_size,shifted_range);

	sprintf(vessel_output,"%s_v_nv%d_nu%d_s%d_vs%d_r%d.dt",output_path.c_str(),num_vessel,num_nuclei_per_vessel,(int)shrink, voxel_size,shifted_range);
	os = new std::ofstream(nuclei_output, std::ios::out | std::ios::binary);
	os2 = new std::ofstream(nuclei_output2, std::ios::out | std::ios::binary);

	// generate some job for worker to process
	int x_dim = (int)pow((float)num_vessel, 1.0/3);
	int y_dim = x_dim;
	int z_dim = num_vessel/(x_dim*y_dim);
	num_nuclei_per_vessel = (num_nuclei_per_vessel*num_vessel)/(x_dim*y_dim*z_dim);
	load_prototype(nuclei_pt.c_str(), vessel_pt.c_str());
	logt("load prototype files", start);
	vector<tuple<float, float, float>> vessel_shifts;
	for(int i=0;i<x_dim;i++){
		for(int j=0;j<y_dim;j++){
			for(int k=0;k<z_dim;k++){
				jobs.push(std::make_tuple(i*vessel_box.high[0], j*vessel_box.high[1], k*vessel_box.high[2]));
				vessel_shifts.push_back(std::make_tuple(i*vessel_box.high[0], j*vessel_box.high[1], k*vessel_box.high[2]));
			}
		}
	}

	pthread_t threads[num_threads];
	for(int i=0;i<num_threads;i++){
		pthread_create(&threads[i], NULL, generate_unit, NULL);
	}
	generate_vessel(vessel_output, vessel_shifts);
	vessel_shifts.clear();
	for(int i = 0; i < num_threads; i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}
	os->close();
	os2->close();
	logt("%ld nucleis are generated for %d vessels", start, global_generated, x_dim*y_dim*z_dim);
	delete os;
	delete vessel;
	for(Polyhedron *n:nucleis){
		delete n;
	}
	nucleis.clear();
}

