/*
 * partition.h
 *
 *  Created on: Nov 18, 2019
 *      Author: teng
 */

#ifndef PARTITION_H_
#define PARTITION_H_
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

// partition the space
void get_mbbs(std::vector<std::string> &input_folders,
		std::vector<aab *> &mbbs, int num_threads);

void persist_tile(std::vector<aab> &tiles, const char *space_path);
void load_space(std::vector<aab> &tiles, const char *space_path);
void partition_data(std::vector<aab> &tiles, const char *output_folder);

/*
 * partition with OCTree
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

/*
 * partition with sorting
 *
 * */

typedef enum partition_direction{
	ROOT,
	ON_X,
	ON_Y,
	ON_Z
} partition_direction;

// sorting partition node
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
	void genTiles(vector<aab> &tiles);
	aab box;
	std::vector<SPNode *> children;
	SPNode *parent;
};

SPNode *build_sort_partition(std::vector<aab*> &mbbs, int num_tiles);


}
#endif /* PARTITION_H_ */
