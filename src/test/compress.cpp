/*****************************************************************************
* Copyright (C) 2011 Adrien Maglo and Clément Courbet
*
* This file is part of PPMC.
*
* PPMC is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* PPMC is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with PPMC.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/

#include <unistd.h>
#include <sys/time.h>
#include <queue>

#include "../storage/tile.h"
#include "../util/util.h"
#include "../spatial/spatial.h"

using namespace CGAL;
using namespace std;
using namespace hispeed;

#define QUEUE_SIZE 100
#define MAX_THREAD_NUM 100

// some shared parameters
string processing_line[MAX_THREAD_NUM];
bool is_working[MAX_THREAD_NUM];
pthread_mutex_t line_lock;
pthread_mutex_t output_lock;
bool stop = false;

std::ofstream *os;
inline void flush_mesh_buffer(vector<HiMesh *> &mesh_buffer, vector<vector<Voxel *>> &voxels){
	assert(mesh_buffer.size()==voxels.size());
	pthread_mutex_lock(&output_lock);
	for(int i=0;i<mesh_buffer.size();i++){
		os->write((char *)&mesh_buffer[i]->dataOffset, sizeof(size_t));
		os->write(mesh_buffer[i]->p_data, mesh_buffer[i]->dataOffset);
		size_t size = voxels[i].size();
		os->write((char *)&size, sizeof(size_t));
		for(Voxel *v:voxels[i]){
			os->write((char *)v->box.min, 3*sizeof(float));
			os->write((char *)v->box.max, 3*sizeof(float));
			os->write((char *)v->core, 3*sizeof(float));
			delete v;
		}
		voxels[i].clear();
		delete mesh_buffer[i];
	}
	pthread_mutex_unlock(&output_lock);
	mesh_buffer.clear();
	voxels.clear();
}

void *mycompress(void *args){
	int id = *(int *)args;
	pthread_mutex_lock(&output_lock);
	log("thread %d is started", id);
	pthread_mutex_unlock(&output_lock);

	// output binary file for the compressed data
	int offset = 0;
	vector<HiMesh *> mesh_buffer;
	vector<vector<Voxel *>> voxels;

	while (!stop||is_working[id]) {
		if(!is_working[id]){
			usleep(10);
			continue;
		}
		MyMesh *mesh = hispeed::get_mesh(processing_line[id], true);
		HiMesh *himesh = new HiMesh(mesh->p_data, mesh->dataOffset);
		himesh->advance_to(100);
		voxels.push_back(himesh->generate_voxels(500));
		delete mesh;
		mesh_buffer.push_back(himesh);
		// if the buffer is full, write the compressed data into binary file
		if(mesh_buffer.size()>=QUEUE_SIZE){
			flush_mesh_buffer(mesh_buffer, voxels);
		}
		processing_line[id].clear();
		is_working[id] = false;
	} // end of while
	// writing last chunk
	flush_mesh_buffer(mesh_buffer, voxels);

	pthread_exit(NULL);
	return NULL;
}


int main(int argc, char** argv) {
	if(argc<2){
		cerr<<"usage: compress input_path output_path [maximum_line]"<<endl;
		exit(0);
	}
	const char *input_path = argv[1];
	const char *output_path = argv[2];
	cerr<<"processing "<<input_path<<" into "<<output_path<<endl;
	int num_threads = hispeed::get_num_threads();
	if(argc>4){
		num_threads = atoi(argv[4]);
	}
	long maximum_lines = LONG_MAX;
	if(argc>3){
		maximum_lines = atoi(argv[3]);
	}
	pthread_t threads[num_threads];
	int id[num_threads];
	for(int i=0;i<num_threads;i++){
		id[i] = i;
		processing_line[i].clear();
		is_working[i] = false;
	}
	// write the compressed data
	os = new std::ofstream(output_path, std::ios::out | std::ios::binary);
	assert(os);
	ifstream is(input_path);
	assert(is);
	long total_filesize = hispeed::file_size(input_path);
	assert(total_filesize>0);

	for(int i=0;i<num_threads;i++){
		pthread_create(&threads[i], NULL, mycompress, (void *)&id[i]);
	}
	struct timeval start_time = get_cur_time();
	std::string input_line;
	int num_objects = 0;
	int next_report = 0;
	long processed_size = 0;

	while (getline(is, input_line)) {
		while(true){
			bool assigned = false;
			for(int i=0;i<num_threads;i++){
				if(is_working[i]==false){
					processing_line[i] = input_line;
					is_working[i] = true;
					assigned = true;
					break;
				}
			}
			if(assigned){
				break;
			}
			usleep(10);
		}
		processed_size += input_line.size()+1;
		num_objects++;
		if(processed_size*100/total_filesize==next_report){
			log("processed %d objects\t(%d\%)", num_objects, next_report++);
		}
		if(num_objects==maximum_lines){
			break;
		}
	}
	stop = true;

	for(int i = 0; i < num_threads; i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}

	os->flush();
	os->close();
	is.close();

	logt("processed %d objects", start_time, num_objects);
	log("total size of compressed data is %ld", hispeed::file_size(output_path));

	pthread_exit(NULL);
	return true;
}

