#include "index.h"

using namespace std;

namespace hispeed{


OctreeNode::OctreeNode(aab b, int _level, long tsize) {
	node_voxel.box = b;
	node_voxel.size = 0;
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
bool OctreeNode::intersects(weighted_aab *object) {
	return node_voxel.box.intersect(object->box);
}

bool OctreeNode::addObject(weighted_aab *object) {
	node_voxel.size += object->size;
	// newly added node
	if(node_voxel.size == object->size){
		// newly added node must be a leaf
		assert(isLeaf);
		objectList.push_back(object);
		return true;
	}
	if (isLeaf) {
		objectList.push_back(object);
		/* Temporary variables */
		if (node_voxel.size > tile_size && canBeSplit) {
			/* Update the center */
			float mid[3];
			float *low = node_voxel.box.min;
			float *high = node_voxel.box.max;
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
			for(weighted_aab *v:objectList) {
				bool inserted = false;
				for (int i = 0; i < 8; i++) {
					if (children[i]->intersects(v)) {
						children[i]->addObject(v);
						totalChildrenSize += v->size;
						inserted = true;
					}
				}
				assert(inserted && "one object must be assigned to at least one child");
			}

			// split the node won't do any help thus we undo
			// this round of split and wait longer
			if (totalChildrenSize >= 2 * (node_voxel.size - 1)) {
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
		} else if (node_voxel.size > 1.5 * tile_size) {
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
			tiles.push_back(node_voxel.box);
		}
	}else{
		for(OctreeNode *n:children){
			n->genTiles(tiles);
		}
	}
}

inline bool update_distance_list(range &d, vector<pair<int, range>> &results){
	bool keep = true;
	int list_size = results.size();
	for(int i = 0;i<list_size;i++){
		if(results[i].second<d){
			keep = false;
		}else if(results[i].second>d){
			results.erase(results.begin()+i);
			list_size--;
			continue;
		}
		i++;
	}
	return keep;
}

void OctreeNode::query_distance(weighted_aab *box, vector<pair<int, range>> &results,
		float min_farthest){
	if(node_voxel.intersect(*box)){
		if(isLeaf){
			for(weighted_aab *obj:objectList){
				if(obj==box){// avoid self comparing
					continue;
				}
				range dis = obj->distance(*box);
				if(dis.closest>min_farthest){
					continue;
				}
				if(update_distance_list(dis, results)){
					results.push_back(pair<int, range>(obj->id, dis));
					if(min_farthest>dis.farthest){
						min_farthest=dis.farthest;
					}
				}
			}
		}else{
			for(OctreeNode *c:this->children){
				c->query_distance(box, results);
			}
		}
	}
}

void OctreeNode::query_intersect(weighted_aab *box, vector<int> &results){
	if(this->node_voxel.intersect(*box)){
		if(this->isLeaf){
			for(weighted_aab *obj:objectList){
				if(obj==box){// avoid self comparing
					continue;
				}
				if(obj->intersect(*box)){
					results.push_back(obj->id);
				}
			}
		}else{
			for(OctreeNode *c:this->children){
				c->query_intersect(box, results);
			}
		}
	}
}

OctreeNode *build_octree(std::vector<weighted_aab*> &voxels, int leaf_size){
	// the main thread build the OCTree with the Minimum Boundary Box
	// get from the data
	aab root_node;
	for(weighted_aab *v:voxels){
		root_node.update(v->box);
	}
	OctreeNode *octree = new OctreeNode(root_node, 0, leaf_size);
	for(weighted_aab *v:voxels){
		octree->addObject(v);
	}
	return octree;
}

}
