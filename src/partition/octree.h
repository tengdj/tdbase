#include <iostream>
#include <sstream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <map>
#include <vector>
#include <cstdlib> 
#include <algorithm>

#include <boost/program_options.hpp>
#include "../geometry/aab.h"
using namespace std;

static int GLOBAL_MAX_LEVEL = 10000000;

namespace hispeed{

extern long tile_size;

class OctreeNode {
	aab box;
public:
	int level;
	bool isLeaf;
	bool canBeSplit;
	OctreeNode* children[8];
	vector<aab*> objectList;

	OctreeNode(aab b, int level);
	~OctreeNode();
	bool addObject(aab *object);
	bool intersects(aab *object);
	void genTiles(vector<aab> &tiles);
};

}
