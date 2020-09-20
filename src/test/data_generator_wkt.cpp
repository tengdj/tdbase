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
Polyhedron vessel;
aab nuclei_box;
aab vessel_box;

float *possibility;

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

void load_prototype(const char *nuclei_path, const char *vessel_path, float shrink, int num_nuclei_per_vessel){
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

	nfile.close();
	vfile.close();
	// update the distance
	int nuclei_num[3];
	for(int i=0;i<3;i++){
		nuclei_num[i] = (int)(vessel_box.max[i]/nuclei_box.max[i]);
	}
	int total_nuclei_num = nuclei_num[0]*nuclei_num[1]*nuclei_num[2];
	float pos_adjust = (float)num_nuclei_per_vessel/total_nuclei_num*1.5;

	possibility = new float[total_nuclei_num];
	float max_distance = DBL_MIN;
	int zero_num = 0;
	for(int i=0;i<nuclei_num[0];i++){
		for(int j=0;j<nuclei_num[1];j++){
			for(int k=0;k<nuclei_num[2];k++){
				float min_dist = DBL_MAX;
				aab cur_box;
				cur_box.min[0] = nuclei_box.max[0]*i;
				cur_box.max[0] = nuclei_box.max[0]*(i+1);
				cur_box.min[1] = nuclei_box.max[1]*j;
				cur_box.max[1] = nuclei_box.max[1]*(j+1);
				cur_box.min[2] = nuclei_box.max[2]*k;
				cur_box.max[2] = nuclei_box.max[2]*(k+1);
				for(Voxel *v:vessel_voxels){
					range r = v->box.distance(cur_box);
					if(min_dist>r.mindist){
						min_dist = r.mindist;
						if(min_dist == 0){
							zero_num++;
							break;
						}
					}
				}
				possibility[i*nuclei_num[1]*nuclei_num[2]+j*nuclei_num[2]+k] = min_dist;
				if(max_distance<min_dist){
					max_distance = min_dist;
				}
			}
		}
	}
	// update the possibility for each voxel
	for(int i=0;i<nuclei_num[0];i++){
		for(int j=0;j<nuclei_num[1];j++){
			for(int k=0;k<nuclei_num[2];k++){
				float cur_possibility = possibility[i*nuclei_num[1]*nuclei_num[2]+j*nuclei_num[2]+k];
				if(cur_possibility!=0){
					// here we control the distribution of the nucleis around vessels
					cur_possibility = (1-cur_possibility/max_distance)*pos_adjust;
					possibility[i*nuclei_num[1]*nuclei_num[2]+j*nuclei_num[2]+k] = cur_possibility;
				}
			}
		}
	}
	for(Voxel *v:vessel_voxels){
		delete v;
	}
	vessel_voxels.clear();
}

Polyhedron *shift_polyhedron(float shift[3], Polyhedron &poly_o){
	Polyhedron *poly = new Polyhedron();
	stringstream ss;
	ss.precision(10);
	ss << poly_o;
	ss >> *poly;
	for(Polyhedron::Vertex_iterator vi=poly->vertices_begin();vi!=poly->vertices_end();vi++){
		Point p = vi->point();
		vi->point() = Point((p[0]+shift[0]), (p[1]+shift[1]), p[2]+shift[2]);
	}
	return poly;
}

/*
 * generate the binary data of nucleis around vessel with
 * a given shift base
 *
 * */
inline int generate_nuclei(float base[3], vector<string> &data1, vector<string> &data2){
	int nuclei_num[3];
	for(int i=0;i<3;i++){
		nuclei_num[i] = (int)(vessel_box.max[i]/nuclei_box.max[i]);
	}
	float shift[3];
	int generated = 0;
	for(int i=0;i<nuclei_num[0];i++){
		for(int j=0;j<nuclei_num[1];j++){
			for(int k=0;k<nuclei_num[2];k++){
				float cur_possibility = possibility[i*nuclei_num[1]*nuclei_num[2]+j*nuclei_num[2]+k];
				if(!hispeed::get_rand_sample(cur_possibility)){
					continue;
				}
				shift[0] = i*nuclei_box.max[0]+base[0];
				shift[1] = j*nuclei_box.max[1]+base[1];
				shift[2] = k*nuclei_box.max[2]+base[2];
				int polyid = hispeed::get_rand_number(nucleis.size()-1);
				int polyid2 = hispeed::get_rand_number(nucleis.size()-1);

				Polyhedron *shifted1 = shift_polyhedron(shift, nucleis[polyid]);
				Polyhedron *shifted2 = shift_polyhedron(shift, nucleis[polyid2]);
				data1.push_back(hispeed::polyhedron_to_wkt(shifted1));
				data2.push_back(hispeed::polyhedron_to_wkt(shifted2));
				delete shifted1;
				delete shifted2;
				generated++;
			}
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
	vector<string> nuclei1;
	vector<string> nuclei2;
	while(!jobs.empty()){
		pthread_mutex_lock(&mylock);

		tuple<float, float, float> job = jobs.front();
		jobs.pop();
		log("%ld jobs left", jobs.size());

		pthread_mutex_unlock(&mylock);
		float base[3] = {get<0>(job),get<1>(job),get<2>(job)};
		int generated = generate_nuclei(base, nuclei1, nuclei2);

		pthread_mutex_lock(&mylock);
		for(string s:nuclei1){
			*os<<s;
		}
		for(string s:nuclei2){
			*os2<<s;
		}
		global_generated += generated;
		pthread_mutex_unlock(&mylock);
		nuclei1.clear();
		nuclei2.clear();
	}
	return NULL;
}

void generate_vessel(const char *path, vector<tuple<float, float, float>> &vessel_shifts){

	ofstream v_os(path);
	for(tuple<float, float, float> tp:vessel_shifts){
		float shift[3] = {get<0>(tp),get<1>(tp),get<2>(tp)};
		Polyhedron *shifted = shift_polyhedron(shift, vessel);
		string s = hispeed::polyhedron_to_wkt(shifted);
		v_os<<s;
	}
	v_os.close();
}


int main(int argc, char **argv){
	string nuclei_pt;
	string vessel_pt;
	string output_path;
	int num_threads = hispeed::get_num_threads();
	int num_nuclei_per_vessel = 10000;
	int num_vessel = 1;
	float shrink = 20;
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
	sprintf(nuclei_output,"%s_n_nv%d_nu%d_s%d.wkt",output_path.c_str(),num_vessel,num_nuclei_per_vessel,(int)shrink);
	sprintf(nuclei_output2,"%s_n2_nv%d_nu%d_s%d.wkt",output_path.c_str(),num_vessel,num_nuclei_per_vessel,(int)shrink);

	sprintf(vessel_output,"%s_v_nv%d_nu%d_s%d.wkt",output_path.c_str(),num_vessel,num_nuclei_per_vessel,(int)shrink);
	os = new std::ofstream(nuclei_output);
	os2 = new std::ofstream(nuclei_output2);

	// generate some job for worker to process
	int x_dim = (int)pow((float)num_vessel, 1.0/3);
	int y_dim = x_dim;
	int z_dim = num_vessel/(x_dim*y_dim);
	num_nuclei_per_vessel = (num_nuclei_per_vessel*num_vessel)/(x_dim*y_dim*z_dim);
	load_prototype(nuclei_pt.c_str(), vessel_pt.c_str(), shrink, num_nuclei_per_vessel);
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
	logt("generated %d vessels",start,vessel_shifts.size());
	vessel_shifts.clear();

	for(int i = 0; i < num_threads; i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}
	os->close();
	os2->close();
	logt("%ld nucleis are generated for %d vessels", start, global_generated, x_dim*y_dim*z_dim);
	delete os;
	delete os2;
}

