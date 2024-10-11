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

#include "aab.h"
#include "RTree.h"
#include "util.h"

#define FillFactor 0.9
#define IndexCapacity 10
#define LeafCapacity 50

using namespace std;
namespace tdbase{


/*
 * OCTree
 * */
class OctreeNode: public weighted_aab{
public:
	long tile_size;
	int level;
	bool isLeaf;
	bool canBeSplit;
	OctreeNode* children[8];
	vector<weighted_aab*> objectList;
	bool isroot(){
		return level==0;
	}
	OctreeNode(aab b, int level, long tsize);
	~OctreeNode();
	bool addObject(weighted_aab *object);

	bool intersects(weighted_aab *object);
	void query_knn(weighted_aab *box, vector<pair<int, range>> &results, float &max_maxdist, const int k=1);
	void query_within(weighted_aab *box, vector<pair<int, range>> &results, const float min_farthest);
	void query_intersect(weighted_aab *box, vector<int> &results);

};
OctreeNode *build_octree(std::vector<weighted_aab*> &mbbs, int num_tiles);

// sorting tree
class SPNode{
	weighted_aab node_voxel;
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
		node_voxel.update(child->node_voxel);
		node_voxel.size += child->node_voxel.size;
		children.push_back(child);
		child->parent = this;
	}
	void add_object(weighted_aab *obj){
		node_voxel.size += obj->size;
		node_voxel.update(*obj);
	}
	bool load(const char *path);
	bool persist(const char *path);

};

SPNode *build_sort_partition(std::vector<weighted_aab*> &mbbs, int num_tiles);

}
#endif /* HISPEED_INDEX_H_ */
