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
#include <spatialindex/SpatialIndex.h>

#define FillFactor 0.9
#define IndexCapacity 10
#define LeafCapacity 50

using namespace std;
namespace hispeed{


/*
 * OCTree
 * */
class OctreeNode {
	weighted_aab node_voxel;
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
	void genTiles(vector<aab> &tiles);
	void query_distance(weighted_aab *box, vector<pair<int, range>> &results, float min_farthest=DBL_MAX);
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
		node_voxel.box.update(child->node_voxel.box);
		node_voxel.size += child->node_voxel.size;
		children.push_back(child);
		child->parent = this;
	}
	void add_object(weighted_aab *obj){
		node_voxel.size += obj->size;
		node_voxel.box.update(obj->box);
	}
	bool load(const char *path);
	bool persist(const char *path);
	void genTiles(vector<aab> &tiles);

};

SPNode *build_sort_partition(std::vector<weighted_aab*> &mbbs, int num_tiles);


class aab_d{
public:
	double low[3];
	double high[3];
	aab_d(aab &b){
		for(int i=0;i<3;i++){
			low[i] = b.min[i];
			high[i] = b.max[i];
		}
	}
};
class MyVisitor : public SpatialIndex::IVisitor{
public:
	std::vector<SpatialIndex::id_type> matches; // contains ids of matching objects
public:
	MyVisitor() {}

	~MyVisitor() {
		matches.clear();
	}

	void visitNode(const SpatialIndex::INode& n) {}
	void visitData(std::string &s) {}

	void visitData(const SpatialIndex::IData& d)
	{
		matches.push_back(d.getIdentifier());
	}

	void visitData(std::vector<const SpatialIndex::IData*>& v) {}
	void visitData(std::vector<uint32_t>& v){}
};

/* Customized data stream to read from a mapping of int to geomery */
class CustomDataStream : public SpatialIndex::IDataStream{
public:
	SpatialIndex::RTree::Data* m_pNext;
	std::vector<aab_d> *shapes;
	std::vector<aab_d>::iterator iter;
	int len;
	SpatialIndex::id_type m_id;
public:
	CustomDataStream(std::vector<aab_d> *inputdata);
	virtual ~CustomDataStream(){
		if (m_pNext != 0) delete m_pNext;
	}

	virtual SpatialIndex::IData* getNext(){
		if (m_pNext == 0) return 0;
		SpatialIndex::RTree::Data* ret = m_pNext;
		m_pNext = 0;
		readNextEntry();
		return ret;
	}

	virtual bool hasNext(){
		return (m_pNext != 0);
	}

	virtual uint32_t size(){
		return len;
	}

	virtual void rewind(){
		if (m_pNext != 0){
			delete m_pNext;
			m_pNext = 0;
		}
		m_id  = 0;
		iter = shapes->begin();
		readNextEntry();
	}
	void readNextEntry();
};

extern bool build_index_geoms(std::vector<aab_d> &geom_mbbs,
							  SpatialIndex::ISpatialIndex* &spidx,
							  SpatialIndex::IStorageManager* &storage);

}
#endif /* HISPEED_INDEX_H_ */
