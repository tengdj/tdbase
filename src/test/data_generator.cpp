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
#include "../spatial/himesh.h"

using namespace hispeed;
using namespace std;

namespace po = boost::program_options;

// 50M for each vessel and the nucleis around it
const int buffer_size = 50*(1<<20);
vector<Polyhedron> nucleis;
vector<vector<Voxel *>> nucleis_voxels;
bool *vessel_taken;
Polyhedron vessel;
aab nuclei_box;
aab vessel_box;


int num_nuclei_per_vessel = 10000;
int num_vessel = 1;
float shrink = 20;
int voxel_size = 400000;

HiMesh *poly_to_himesh(Polyhedron &poly){
	stringstream ss;
	ss<<poly;
	MyMesh *mesh = hispeed::get_mesh(ss.str(), true);
	HiMesh *himesh = new HiMesh(mesh->p_data, mesh->dataOffset);
	himesh->advance_to(100);
	delete mesh;
	return himesh;
}

MyMesh *poly_to_mesh(Polyhedron poly){
	stringstream ss;
	ss<<poly;
	return hispeed::get_mesh(ss.str(), true);
}

void load_prototype(const char *nuclei_path, const char *vessel_path){
	std::ifstream nfile(nuclei_path);
	std::ifstream vfile(vessel_path);
	string input_line;
	vector<Voxel *> vessel_voxels;
	if(std::getline(vfile, input_line)){
		hispeed::replace_bar(input_line);
		stringstream ss;
		ss<<input_line;
		ss >> vessel;
		aab tmpb;
		for(Polyhedron::Vertex_iterator vi=vessel.vertices_begin();vi!=vessel.vertices_end();vi++){
			Point p = vi->point();
			tmpb.update(p[0], p[1], p[2]);
		}
		for(Polyhedron::Vertex_iterator vi=vessel.vertices_begin();vi!=vessel.vertices_end();vi++){
			Point p = vi->point();
			vi->point() = Point(p[0]-tmpb.min[0], p[1]-tmpb.min[1], p[2]-tmpb.min[2]);
		}
		tmpb.max[0] -= tmpb.min[0];
		tmpb.max[1] -= tmpb.min[1];
		tmpb.max[2] -= tmpb.min[2];
		tmpb.min[0] = 0;
		tmpb.min[1] = 0;
		tmpb.min[2] = 0;
		vessel_box.update(tmpb);
		HiMesh *himesh = poly_to_himesh(vessel);
		vessel_voxels = himesh->generate_voxels(100);

		delete himesh;
	}else{
		assert(false&&"error reading the vessel file");
	}

	while(std::getline(nfile, input_line)){
		hispeed::replace_bar(input_line);
		stringstream ss;
		ss<<input_line;
		Polyhedron poly;
		ss >> poly;
		aab tmpb;
		for(Polyhedron::Vertex_iterator vi=poly.vertices_begin();vi!=poly.vertices_end();vi++){
			Point p = vi->point();
			Point np(p[0]/shrink, p[1]/shrink, p[2]/shrink);
			vi->point() = np;
			tmpb.update(np[0], np[1], np[2]);
		}
		for(Polyhedron::Vertex_iterator vi=poly.vertices_begin();vi!=poly.vertices_end();vi++){
			Point p = vi->point();
			vi->point() = Point((p[0]-tmpb.min[0]), (p[1]-tmpb.min[1]), p[2]-tmpb.min[2]);
		}
		tmpb.max[0] -= tmpb.min[0];
		tmpb.max[1] -= tmpb.min[1];
		tmpb.max[2] -= tmpb.min[2];
		tmpb.min[2] = 0;
		tmpb.min[1] = 0;
		tmpb.min[2] = 0;
		nuclei_box.update(tmpb);
		nucleis.push_back(poly);
		HiMesh *himesh = poly_to_himesh(poly);
		nucleis_voxels.push_back(himesh->generate_voxels(500));
		delete himesh;
	}


	int nuclei_num[3];
	for(int i=0;i<3;i++){
		nuclei_num[i] = (int)(vessel_box.max[i]/nuclei_box.max[i]);
	}

	int total_slots = nuclei_num[0]*nuclei_num[1]*nuclei_num[2];
	vessel_taken = new bool[total_slots];
	for(Voxel *v:vessel_voxels){
		int xstart = v->box.min[0]/nuclei_num[0];
		int xend = v->box.max[0]/nuclei_num[0];
		int ystart = v->box.min[1]/nuclei_num[1];
		int yend = v->box.max[1]/nuclei_num[1];
		int zstart = v->box.min[2]/nuclei_num[2];
		int zend = v->box.max[2]/nuclei_num[2];
		for(int z=zstart;z<=zend;z++){
			for(int y=ystart;y<=yend;y++){
				for(int x=xstart;x<=xend;x++){
					vessel_taken[z*nuclei_num[0]*nuclei_num[1]+y*nuclei_num[0]+x] = true;
				}
			}
		}

	}

	nfile.close();
	vfile.close();

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
			box_tmp[i] = v->box.min[i]+shift[i];
		}
		memcpy(data+offset, (char *)box_tmp, 3*sizeof(float));
		offset += 3*sizeof(float);

		for(int i=0;i<3;i++){
			box_tmp[i] = v->box.max[i]+shift[i];
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
		nuclei_num[i] = (int)(vessel_box.max[i]/nuclei_box.max[i]);
	}
	float shift[3];
	int generated = 0;
	int total_slots = nuclei_num[0]*nuclei_num[1]*nuclei_num[2];
	bool *taken = new bool[total_slots];

	assert(total_slots>num_nuclei_per_vessel);

	while(++generated<num_nuclei_per_vessel){
		int idx = hispeed::get_rand_number(total_slots);
		while(taken[idx]||vessel_taken[idx]){
			idx = hispeed::get_rand_number(total_slots);
		}
		taken[idx] = true;

		int z = idx/(nuclei_num[0]*nuclei_num[1]);
		int y = (idx%(nuclei_num[0]*nuclei_num[1]))/nuclei_num[0];
		int x = (idx%(nuclei_num[0]*nuclei_num[1]))%nuclei_num[0];

		shift[0] = x*nuclei_box.max[0]+base[0];
		shift[1] = y*nuclei_box.max[1]+base[1];
		shift[2] = z*nuclei_box.max[2]+base[2];

		int polyid = hispeed::get_rand_number(nucleis.size()-1);
		organize_data(nucleis[polyid], nucleis_voxels[polyid], shift, data, offset);
		{
			float shift2[3];
			shift2[0] = shift[0]+nuclei_box.max[0]*(hispeed::get_rand_number(100)*1.0)/100.0*(hispeed::get_rand_sample(50)?1:-1);
			shift2[1] = shift[1]+nuclei_box.max[1]*(hispeed::get_rand_number(100)*1.0)/100.0*(hispeed::get_rand_sample(50)?1:-1);
			shift2[2] = shift[2]+nuclei_box.max[2]*(hispeed::get_rand_number(100)*1.0)/100.0*(hispeed::get_rand_sample(50)?1:-1);

			int polyid2 = hispeed::get_rand_number(nucleis.size()-1);
			organize_data(nucleis[polyid2], nucleis_voxels[polyid2], shift2, data2, offset2);
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
	while(!jobs.empty()){
		pthread_mutex_lock(&mylock);

		tuple<float, float, float> job = jobs.front();
		jobs.pop();
		log("%ld jobs left", jobs.size());

		pthread_mutex_unlock(&mylock);
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
	char *data = new char[vessel_shifts.size()*100000*2];
	size_t offset = 0;
	HiMesh *himesh = poly_to_himesh(vessel);
	vector<Voxel *> voxels = himesh->generate_voxels(voxel_size);
	for(tuple<float, float, float> tp:vessel_shifts){
		float shift[3] = {get<0>(tp),get<1>(tp),get<2>(tp)};
		organize_data(vessel, voxels, shift, data, offset);
	}
	ofstream *v_os = new std::ofstream(path, std::ios::out | std::ios::binary);
	v_os->write(data, offset);
	v_os->close();
	delete v_os;
	log("%d voxels are generated for the vessel",voxels.size());
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
	sprintf(nuclei_output,"%s_n_nv%d_nu%d_s%d_vs%d.dt",output_path.c_str(),num_vessel,num_nuclei_per_vessel,(int)shrink, voxel_size);
	sprintf(nuclei_output2,"%s_n2_nv%d_nu%d_s%d_vs%d.dt",output_path.c_str(),num_vessel,num_nuclei_per_vessel,(int)shrink, voxel_size);

	sprintf(vessel_output,"%s_v_nv%d_nu%d_s%d_vs%d.dt",output_path.c_str(),num_vessel,num_nuclei_per_vessel,(int)shrink, voxel_size);
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
				jobs.push(std::make_tuple(i*vessel_box.max[0], j*vessel_box.max[1], k*vessel_box.max[2]));
				vessel_shifts.push_back(std::make_tuple(i*vessel_box.max[0], j*vessel_box.max[1], k*vessel_box.max[2]));
			}
		}
	}
	if(num_threads>jobs.size()){
		num_threads = jobs.size();
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
}

