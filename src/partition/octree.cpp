#include "partition.h"

using namespace std;

namespace hispeed{


OctreeNode::OctreeNode(aab b, int _level, long tsize) {
	box = b;
	level = _level;
	isLeaf = true;
	canBeSplit = true;
	tile_size = tsize;
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
			children[0] = new OctreeNode(aab(low[0], low[1], low[2], mid[0], mid[1], mid[2]), level + 1, tile_size);
			children[1] = new OctreeNode(aab(mid[0], low[1], low[2], high[0], mid[1], mid[2]), level + 1, tile_size);
			children[2] = new OctreeNode(aab(low[0], mid[1], low[2], mid[0], high[1], mid[2]), level + 1, tile_size);
			children[3] = new OctreeNode(aab(mid[0], mid[1], low[2], high[0], high[1], mid[2]), level + 1, tile_size);
			children[4] = new OctreeNode(aab(low[0], low[1], mid[2], mid[0], mid[1], high[2]), level + 1, tile_size);
			children[5] = new OctreeNode(aab(mid[0], low[1], mid[2], high[0], mid[1], high[2]), level + 1, tile_size);
			children[6] = new OctreeNode(aab(low[0], mid[1], mid[2], mid[0], high[1], high[2]), level + 1, tile_size);
			children[7] = new OctreeNode(aab(mid[0], mid[1], mid[2], high[0], high[1], high[2]), level + 1, tile_size);

			// assign objects to children
			int totalChildrenSize = 0;
			for(aab *b:objectList) {
				bool inserted = false;
				for (int i = 0; i < 8; i++) {
					if (children[i]->intersects(b)) {
						children[i]->addObject(b);
						totalChildrenSize += b->weight;
						inserted = true;
					}
				}
				assert(inserted && "one object must be assigned to at least one child");
			}

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
			tiles.push_back(box);
		}
	}else{
		for(OctreeNode *n:children){
			n->genTiles(tiles);
		}
	}
}

OctreeNode *build_ctree(std::vector<aab*> &mbbs, int num_tiles){
	// the main thread build the OCTree with the Minimum Boundary Box
	// get from the data
	aab root_node;
	long total_size = 0;
	for(aab *b:mbbs){
		root_node.update(*b);
		total_size += b->weight;
	}
	OctreeNode *octree = new OctreeNode(root_node, 0, total_size/num_tiles);
	for(aab *b:mbbs){
		octree->addObject(b);
	}
	return octree;
}

}