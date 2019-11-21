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

#include "../PPMC/ppmc.h"
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
inline void flush_mesh_buffer(vector<MyMesh *> &mesh_buffer){
	pthread_mutex_lock(&output_lock);
	for(MyMesh *mesh:mesh_buffer){
		os->write((char *)&mesh->dataOffset, sizeof(size_t));
		os->write(mesh->p_data, mesh->dataOffset);
		delete mesh;
	}
	pthread_mutex_unlock(&output_lock);
	mesh_buffer.clear();
}

void *compress(void *args){
	int id = *(int *)args;
	pthread_mutex_lock(&output_lock);
	cout<<"thread "<<id<<" is started"<<endl;
	pthread_mutex_unlock(&output_lock);

	// output binary file for the compressed data
	int offset = 0;
	vector<MyMesh *> mesh_buffer;

	while (!stop||is_working[id]) {
		if(!is_working[id]){
			usleep(10);
			continue;
		}
		mesh_buffer.push_back(hispeed::get_mesh(processing_line[id], true));
		// if the buffer is full, write the compressed data into binary file
		if(mesh_buffer.size()>=QUEUE_SIZE){
			flush_mesh_buffer(mesh_buffer);
		}
		processing_line[id].clear();
		is_working[id] = false;
	} // end of while
	// writing last chunk
	flush_mesh_buffer(mesh_buffer);

	pthread_exit(NULL);
	return NULL;
}


int main(int argc, char** argv) {
	if(argc<2){
		cerr<<"usage: compress input_path output_path"<<endl;
		exit(0);
	}
	const char *input_path = argv[1];
	const char *output_path = argv[2];
	cerr<<"processing "<<input_path<<" into "<<output_path<<endl;
	int num_threads = hispeed::get_num_threads();
	assert(num_threads>0 && num_threads<MAX_THREAD_NUM);
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
		int rc = pthread_create(&threads[i], NULL, compress, (void *)&id[i]);
		if (rc) {
			cout << "Error:unable to create thread," << rc << endl;
			exit(-1);
		}
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
			cerr<<"processed "<<num_objects<<" objects\t("<<next_report<<"%)"<<endl;
			next_report++;
		}
	}
	stop = true;

	for(int i = 0; i < num_threads; i++ ){
		void *status;
		int rc = pthread_join(threads[i], &status);
		if (rc) {
			cout << "Error:unable to join," << rc << endl;
			exit(-1);
		}
		cerr << "Main: completed thread id :" << i ;
		cerr << "  exiting with status :" << status << endl;
	}

	os->flush();
	os->close();
	is.close();


    std::cerr <<"processed "<<num_objects<<" objects in "<<get_time_elapsed(start_time)/1000<<" seconds"<<endl;
	std::cerr <<"total size of compressed data is "<< hispeed::file_size(output_path) << std::endl; // the total size

	cerr << "Main: program exiting." << endl;
	pthread_exit(NULL);

	return true;
}

