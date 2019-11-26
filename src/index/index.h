/*
 * index.h
 *
 *  Created on: Nov 18, 2019
 *      Author: teng
 */

#ifndef HISPEED_INDEX_H_
#define HISPEED_INDEX_H_
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <map>
#include <vector>
#include <cstdlib>
#include <algorithm>

#include "../geometry/aab.h"
#include "../PPMC/mymesh.h"
#include "../spatial/spatial.h"

using namespace std;
namespace hispeed{

/*
 * OCTree
 * */
class OctreeNode {
	Voxel node_voxel;
public:
	long tile_size;
	int level;
	bool isLeaf;
	bool canBeSplit;
	OctreeNode* children[8];
	vector<Voxel*> objectList;

	OctreeNode(aab b, int level, long tsize);
	~OctreeNode();
	bool addObject(Voxel *object);
	bool intersects(Voxel *object);
	void genTiles(vector<aab> &tiles);
};
OctreeNode *build_ctree(std::vector<Voxel*> &mbbs, int num_tiles);

// sorting tree
class SPNode{
	Voxel node_voxel;
	std::vector<SPNode *> children;
	SPNode *parent;
public:
	SPNode(){
		parent = NULL;
	}
	~SPNode(){
		for(SPNode *s:children){
			delete s;
		}
		children.clear();
	}
	void add_child(SPNode *child){
		node_voxel.box.update(child->node_voxel.box);
		node_voxel.size += child->node_voxel.size;
		children.push_back(child);
		child->parent = this;
	}
	void add_object(Voxel *obj){
		node_voxel.size += obj->size;
		node_voxel.box.update(obj->box);
	}
	bool load(const char *path);
	bool persist(const char *path);
	void genTiles(vector<aab> &tiles);

};

SPNode *build_sort_partition(std::vector<Voxel*> &mbbs, int num_tiles);


}
#endif /* HISPEED_INDEX_H_ */
