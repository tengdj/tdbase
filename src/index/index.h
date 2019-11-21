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
	aab box;
public:
	long tile_size;
	int level;
	bool isLeaf;
	bool canBeSplit;
	OctreeNode* children[8];
	vector<aab*> objectList;

	OctreeNode(aab b, int level, long tsize);
	~OctreeNode();
	bool addObject(aab *object);
	bool intersects(aab *object);
	void genTiles(vector<aab> &tiles);
};
OctreeNode *build_ctree(std::vector<aab*> &mbbs, int num_tiles);

// sorting tree
class SPNode{
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
	bool load(const char *path);
	bool persist(const char *path);
	void genTiles(vector<aab> &tiles);
	aab box;
	std::vector<SPNode *> children;
	SPNode *parent;
};

SPNode *build_sort_partition(std::vector<aab*> &mbbs, int num_tiles);


}
#endif /* HISPEED_INDEX_H_ */
