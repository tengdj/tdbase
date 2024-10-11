/*
 * sortpartition.cpp
 *
 *  Created on: Nov 19, 2019
 *      Author: teng
 *
 *  partition the space with sorted points
 *
 */

#include "index.h"

#include<iostream>
#include<fstream>
using namespace std;

namespace tdbase{

bool SPNode::load(const char *path){
	ifstream is(path, ios::out | ios::binary);
	if(!is) {
		log("Cannot open file!");
		return false;
	}
	long part_x, part_y, part_z;
	is.read((char *)&part_x, sizeof(long));
	is.read((char *)&part_y, sizeof(long));
	is.read((char *)&part_z, sizeof(long));
	for(int x=0;x<part_x;x++){
		SPNode *cur_x = new SPNode();
		for(int y=0;y<part_y;y++){
			SPNode *cur_y = new SPNode();
			for(int z=0;z<part_z;z++){
				SPNode *cur_z = new SPNode();
				is.read((char *)cur_z->node_voxel.low, sizeof(float)*3);
				is.read((char *)cur_z->node_voxel.high, sizeof(float)*3);
				is.read((char *)&cur_z->node_voxel.size, sizeof(int));
				cur_y->add_child(cur_z);
			}
			cur_x->add_child(cur_y);
		}
		add_child(cur_x);
	}
	is.close();
	return true;
};
bool SPNode::persist(const char *path){
	// only the root node can call the persist function;
	long part_x = children.size();
	if(part_x==0){
		return false;
	}
	long part_y = children[0]->children.size();
	if(part_y==0){
		return false;
	}
	long part_z = children[0]->children[0]->children.size();
	if(part_z==0){
		return false;
	}

	ofstream os(path, ios::out | ios::binary);
	if(!os) {
		log("Cannot open file!");
		return false;
	}
	os.write((char *)&part_x, sizeof(long));
	os.write((char *)&part_y, sizeof(long));
	os.write((char *)&part_z, sizeof(long));
	for(int x=0;x<part_x;x++){
		for(int y=0;y<part_y;y++){
			for(int z=0;z<part_z;z++){
				os.write((char *)children[x]->children[y]->children[z]->node_voxel.low, sizeof(float)*3);
				os.write((char *)children[x]->children[y]->children[z]->node_voxel.high, sizeof(float)*3);
				os.write((char *)&children[x]->children[y]->children[z]->node_voxel.size, sizeof(int));
			}
		}
	}
	os.close();
	return true;
}

bool compareAAB_x(weighted_aab *a1, weighted_aab *a2){
	return a1->low[0]<a2->low[0];
}
bool compareAAB_y(weighted_aab *a1, weighted_aab *a2){
	return a1->low[1]<a2->low[1];
}
bool compareAAB_z(weighted_aab *a1, weighted_aab *a2){
	return a1->low[2]<a2->low[2];
}

SPNode *build_sort_partition(std::vector<weighted_aab*> &voxels, int num_tiles){
	// create the root node
	SPNode *root = new SPNode();

	// get the number of cuts in each dimension
	long part_x = std::max((int)std::pow(num_tiles, 1.0/3),1);
	long part_y = part_x;
	long part_z = std::max((int)((num_tiles/part_x)/part_y),1);

	// sort the mbbs in X dimension
	std::sort(voxels.begin(), voxels.end(), compareAAB_x);
	long length_r = voxels.size(); // length of the root nodes
	long step_x = (length_r+part_x-1)/part_x;
	for(long x=0;x<length_r;x+=step_x){
		long length_x = std::min(length_r, x+step_x) - x;
		SPNode *cur_x = new SPNode();

		// process the slide cut in X dimension
		std::vector<weighted_aab *> tmp_x;
		for(int t=x;t<x+length_x;t++){
			tmp_x.push_back(voxels[t]);
		}
		std::sort(tmp_x.begin(), tmp_x.end(), compareAAB_y);
		long step_y = (length_x+part_y-1)/part_y;
		for(long y=0;y<length_x;y+=step_y){
			long length_y = std::min(length_x, y+step_y) - y;
			SPNode *cur_y = new SPNode();
			// process the slice cut in y dimension
			std::vector<weighted_aab *> tmp_y;
			for(int t=y;t<y+length_y;t++){
				tmp_y.push_back(tmp_x[t]);
			}
			std::sort(tmp_y.begin(), tmp_y.end(), compareAAB_z);
			long step_z = (length_y+part_z-1)/part_z;
			for(long z=0;z<length_y;z+=step_z){
				long length_z = std::min(length_y, z+step_z) - z;
				// generate the leaf node, update its MBB
				// with the objects assigned to it
				SPNode *cur_z = new SPNode();
				for(int t=z;t<z+length_z;t++){
					cur_z->add_object(tmp_y[t]);
				}
				cur_y->add_child(cur_z);
			}
			tmp_y.clear();
			cur_x->add_child(cur_y);
		}
		tmp_x.clear();
		root->add_child(cur_x);
	}
	return root;
}

}
