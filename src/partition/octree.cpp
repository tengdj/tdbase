#include "octree.h"

using namespace std;

namespace hispeed{

// Global variables (used by QuadtreeNode)
long tile_size = 1<<20;

OctreeNode::OctreeNode(aab b, int _level) {
	box = b;
	level = _level;
	isLeaf = true;
	canBeSplit = true;
}

OctreeNode::~OctreeNode() {
	objectList.clear();
	if(!isLeaf){
		for(OctreeNode *c:children){
			delete c;
		}
	}
}

/* Check if the node MBB intersects with the object MBB */
bool OctreeNode::intersects(aab *object) {
	return box.intersect(object);
}

bool OctreeNode::addObject(aab *object) {
	box.weight += object->weight;
	// newly added node
	if(box.weight == object->weight){
		// newly added node must be a leaf
		assert(isLeaf);
		objectList.push_back(object);
		return true;
	}
	if (isLeaf) {
		objectList.push_back(object);
		/* Temporary variables */
		if (box.weight > tile_size && canBeSplit) {
			/* Update the center */
			float mid[3];
			float *low = box.min;
			float *high = box.max;
			for (int i = 0; i < 3; i++) {
				mid[i] = (low[i] + high[i]) / 2;
			}

			// Split the node to 8 nodes equally centered with the middle point
			// no objects are assigned to the children yet
			children[0] = new OctreeNode(aab(low[0], low[1], low[2], mid[0], mid[1], mid[2]), level + 1);
			children[1] = new OctreeNode(aab(mid[0], low[1], low[2], high[0], mid[1], mid[2]), level + 1);
			children[2] = new OctreeNode(aab(low[0], mid[1], low[2], mid[0], high[1], mid[2]), level + 1);
			children[3] = new OctreeNode(aab(mid[0], mid[1], low[2], high[0], high[1], mid[2]), level + 1);
			children[4] = new OctreeNode(aab(low[0], low[1], mid[2], mid[0], mid[1], high[2]), level + 1);
			children[5] = new OctreeNode(aab(mid[0], low[1], mid[2], high[0], mid[1], high[2]), level + 1);
			children[6] = new OctreeNode(aab(low[0], mid[1], mid[2], mid[0], high[1], high[2]), level + 1);
			children[7] = new OctreeNode(aab(mid[0], mid[1], mid[2], high[0], high[1], high[2]), level + 1);

			// assign objects to children
			int totalChildrenSize = 0;
			int good = 0, bad = 0;
			for(aab *b:objectList) {
				int added = 0;
				for (int i = 0; i < 8; i++) {
					if (children[i]->intersects(b)) {
						children[i]->addObject(b);
						totalChildrenSize += b->weight;
						added++;
					}
				}
				if(added==0){
					bad++;
					//cerr<<"not belong to any space "<<*b<<endl;
				}else{
					good++;
					//cerr<<added<<endl;
				}
			}
			cerr<<good<<" "<<bad<<endl;


			// split the node won't do any help thus we undo
			// this round of split and wait longer
			if (totalChildrenSize >= 8 * (box.weight - 1)) {
				isLeaf = true;
				canBeSplit = false;
				for (int i = 0; i < 8; i++) {
					delete children[i];
					children[i] = NULL;
				}
			} else {
				objectList.clear();
				isLeaf = false;
			}
		} else if (box.weight > 1.5 * tile_size) {
			canBeSplit = true;
		}

	} else {
		for (int i = 0; i < 8; i++) {
			if (children[i]->intersects(object)) {
				children[i]->addObject(object);
			}
		}
	}
	return true;
}

void OctreeNode::genTiles(vector<aab> &tiles){
	if(isLeaf){
		if(objectList.size()>0){
			cerr<<this->level<<endl;
			tiles.push_back(box);
		}
	}else{
		for(OctreeNode *n:children){
			n->genTiles(tiles);
		}
	}
}



}
