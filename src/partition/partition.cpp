/*
 * partition.cpp
 *
 *  Created on: Nov 18, 2019
 *      Author: teng
 *  partition the input dataset into tiles. Each tile contains
 *  more or less the same number of computational units (edges or facets)
 *
 */
#include<pthread.h>
#include <queue>
#include <boost/algorithm/string/replace.hpp>

#include "partition.h"
#include "../util/util.h"

using namespace std;

namespace hispeed{

const int LOCAL_QUEUE_SIZE=10;

struct unit_args{
	int id;
	std::queue<string> *input_lines;
	std::vector<weighted_aab *> *voxels;
	bool *complete;
};
pthread_mutex_t lock;

void *get_mbbs_unit(void *arg){
	struct unit_args *arg_val = (struct unit_args *)arg;

	log("thread %d is started", arg_val->id);

	vector<string> local_input_lines;
	vector<weighted_aab*> local_voxels;
	int local_processed = 0;
	while(!(*arg_val->complete)||!arg_val->input_lines->empty()){
		// fetch a certain number of lines
		// from the queue
		pthread_mutex_lock(&lock);
		for(int i=0;i<LOCAL_QUEUE_SIZE;i++){
			if(arg_val->input_lines->empty()){
				break;
			}
			local_input_lines.push_back(arg_val->input_lines->front());
			arg_val->input_lines->pop();
		}
		pthread_mutex_unlock(&lock);
		if(local_input_lines.empty()){
			usleep(10);
			continue;
		}
		// process the polyhedrons
		for(string input_line:local_input_lines){
			MyMesh *mesh = get_mesh(input_line);
			if(mesh==NULL){
				continue;
			}
			weighted_aab *voxel = new weighted_aab();
			voxel->box = aab(mesh->bbMin[0], mesh->bbMin[1], mesh->bbMin[2],
					mesh->bbMax[0], mesh->bbMax[1], mesh->bbMax[2]);
			voxel->size = mesh->size_of_halfedges()/2;
			local_voxels.push_back(voxel);
			delete mesh;
		}
		local_input_lines.clear();
	}

	pthread_mutex_lock(&lock);
	arg_val->voxels->insert(arg_val->voxels->end(), local_voxels.begin(), local_voxels.end());
	local_voxels.clear();
	pthread_mutex_unlock(&lock);
}

void get_voxels(std::vector<std::string> &input_folders, std::vector<weighted_aab *> &voxels,
		const int num_threads, const int sample_rate){
	bool complete = false;
	// get file list
	std::vector<string> files;
	for(string f:input_folders){
		hispeed::list_files(f.c_str(), files);
	}
	const long total_filesize = hispeed::file_size(files);
	long processed_filesize = 0;

	std::queue<string> input_lines;

	// start threads
	pthread_t threads[num_threads];
	struct unit_args arg[num_threads];
	for(int i=0;i<num_threads;i++){
		arg[i].id = i;
		arg[i].input_lines = &input_lines;
		arg[i].voxels = &voxels;
		arg[i].complete = &complete;
		pthread_create(&threads[i], NULL, get_mbbs_unit, (void *)&arg[i]);
	}

	// the main thread read from files and store the lines
	// in a buffer queue with a maximum length
	int read_lines = 0;
	int processed_lines = 0;
	vector<string> local_input_lines;
	long next_report = 1;
	for(string path:files){
		std::ifstream is(path.c_str());
		string input_line;
		while(getline(is, input_line)){
			read_lines++;
			processed_filesize += input_line.size()+1;
			// with a certain sample rate
			if(!hispeed::get_rand_sample(sample_rate)){
				continue;
			}
			processed_lines++;
			// here we use a local buffer to avoid
			// frequent lock request
			local_input_lines.push_back(input_line);
			if(local_input_lines.size()==LOCAL_QUEUE_SIZE){
				pthread_mutex_lock(&lock);
				for(string s:local_input_lines){
					input_lines.push(s);
				}
				pthread_mutex_unlock(&lock);
				local_input_lines.clear();
			}
			// wait if the queue length is too long
			while(input_lines.size()>num_threads*LOCAL_QUEUE_SIZE){
				usleep(10);
			}
			if(processed_filesize*100/total_filesize==next_report){
				log("read %d processed %d objects\t(%ld\%)", read_lines, processed_lines, next_report);
				next_report += 2;
			}
		}
	}
	complete = true;
	// wait for threads to stop
	for(int i = 0; i < num_threads; i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}
}
void persist_tile(std::vector<aab> &tiles, const char *prefix){
	for(int i=0;i<tiles.size();i++){
		Polyhedron *poly = make_cube(tiles[i]);
		stringstream ss;
		ss<<prefix<<i<<".off";
		write_polyhedron(poly, ss.str().c_str());
	}
}

void partition_data(std::vector<aab> &tiles, const char *output_folder){
	;
}


}

