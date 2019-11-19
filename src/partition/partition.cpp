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
#include <boost/algorithm/string/replace.hpp>

#include "partition.h"
#include "../util/util.h"
#include <queue>

using namespace std;

namespace hispeed{

std::vector<aab *> objects;
pthread_mutex_t lock;
std::queue<string> input_lines;

const int local_buffer_size=10;
bool complete = false;

int processed = 0;

void *partition_unit(void *arg){
	pthread_mutex_lock(&lock);
	cerr<<"thread "<<*(int*)arg<<" is started"<<endl;
	pthread_mutex_unlock(&lock);

	vector<string> local_input_lines;
	vector<aab*> local_objects;
	int local_processed = 0;
	while(!complete||!input_lines.empty()){
		pthread_mutex_lock(&lock);
		for(int i=0;i<local_buffer_size;i++){
			if(input_lines.empty()){
				break;
			}
			local_input_lines.push_back(input_lines.front());
			input_lines.pop();
		}
		pthread_mutex_unlock(&lock);
		for(string input_line:local_input_lines){
			boost::replace_all(input_line, "|", "\n");
			MyMesh *mesh = get_mesh(input_line);
			if(mesh==NULL){
				continue;
			}
			aab *box = new aab(mesh->bbMin[0], mesh->bbMin[1], mesh->bbMin[2],
					mesh->bbMax[0], mesh->bbMax[1], mesh->bbMax[2]);
			box->weight = mesh->size_of_facets();
			local_objects.push_back(box);
			delete mesh;
		}
		local_input_lines.clear();
	}

	pthread_mutex_lock(&lock);
	objects.insert(objects.end(), local_objects.begin(), local_objects.end());
	local_objects.clear();
	pthread_mutex_unlock(&lock);
}

void partition_space(std::vector<std::string> input_folders,
		std::vector<aab> &tiles, int num_threads, int num_tiles){
	complete = false;
	processed = 0;
	std::vector<string> files;
	for(string f:input_folders){
		hispeed::list_files(f.c_str(), files);
	}

	pthread_t threads[num_threads];
	int id[num_threads];
	for(int i=0;i<num_threads;i++){
		id[i] = i;
	}

	for(int i=0;i<num_threads;i++){
		int rc = pthread_create(&threads[i], NULL, partition_unit, (void *)&id[i]);
		if (rc) {
			cout << "Error:unable to create thread," << rc << endl;
			exit(-1);
		}
	}

	int index = 0;
	vector<string> local_input_lines;
	for(string path:files){
		std::ifstream is(path.c_str());
		string input_line;
		while(getline(is, input_line)){
			local_input_lines.push_back(input_line);
			if(local_input_lines.size()==local_buffer_size){
				pthread_mutex_lock(&lock);
				for(string s:local_input_lines){
					input_lines.push(s);
				}
				pthread_mutex_unlock(&lock);
				local_input_lines.clear();
			}
			if(++index%100==0){
				cerr<<"processed "<<index<<" "<<input_lines.size()<<endl;
			}
			while(input_lines.size()>num_threads*local_buffer_size){
				usleep(10);
			}
		}
	}
	complete = true;

	pthread_mutex_lock(&lock);
	cerr<<"waiting for workers to stop"<<endl;
	pthread_mutex_unlock(&lock);

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

	aab root_node;
	long total_size = 0;
	for(aab *b:objects){
		root_node.update(*b);
		total_size += b->weight;
	}
	cerr<<root_node<<endl;
	cerr<<objects.size()<<" "<<total_size<<endl;
	tile_size = total_size/num_tiles;
	OctreeNode *octree = new OctreeNode(root_node, 0);
	for(aab *b:objects){
		octree->addObject(b);
	}
	octree->genTiles(tiles);
	delete octree;
	for(aab *b:objects){
		delete b;
	}
	objects.clear();
}
void persist_tile(std::vector<aab> &tiles, const char *space_path){
	;
}

void load_space(std::vector<aab> &tiles, const char *space_path){
	;
}
void partition_data(std::vector<aab> &tiles, const char *output_folder){
	;
}


}

