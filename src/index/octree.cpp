#include "index.h"

using namespace std;

namespace hispeed{


OctreeNode::OctreeNode(aab b, int _level, long tsize) {
	set_box(b);
	size = 0;
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
	return intersect(*object);
}

bool OctreeNode::addObject(weighted_aab *object) {
	size += object->size;
	// newly added node
	if(size == object->size){
		// newly added node must be a leaf
		assert(isLeaf);
		objectList.push_back(object);
		return true;
	}
	if (isLeaf) {
		objectList.push_back(object);
		/* Temporary variables */
		if (size > tile_size && canBeSplit) {
			/* Update the center */
			float mid[3];
			float *low = min;
			float *high = max;
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
			if (totalChildrenSize >= 2 * (size - 1)) {
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
		} else if (size > 1.5 * tile_size) {
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

void OctreeNode::query_nn(weighted_aab *box, vector<pair<int, range>> &candidates, float & min_maxdist){
	range dis = distance(*box);

	//current node possibly covers nearest objects dis.mindist<=min_maxdist
	if(dis.mindist<=min_maxdist){
		if(isLeaf){
			for(weighted_aab *obj:objectList){
				if(obj==box){// avoid self comparing
					continue;
				}
				range objdis = obj->distance(*box);
				if(objdis.mindist<=min_maxdist){

					// check each candidate in results list to see if this one
					// can be kept, or some candidates need be removed from list
					if(update_distance_list(objdis, candidates)){
						candidates.push_back(pair<int, range>(obj->id, objdis));
						// update MINMAXDIST if encountered
						if(min_maxdist>objdis.maxdist){
							min_maxdist=objdis.maxdist;
						}
					}
				}

			}
		}else{
			for(OctreeNode *c:this->children){
				c->query_nn(box, candidates, min_maxdist);
			}
		}
	}
}

void OctreeNode::query_within(weighted_aab *box, vector<pair<int, range>> &results, const float max_dist){
	range dis = distance(*box);

	if(dis.mindist<=max_dist){
		if(isLeaf){
			for(weighted_aab *obj:objectList){
				if(obj==box){// avoid self comparing
					continue;
				}
				range objdis = obj->distance(*box);
				if(objdis.mindist<=max_dist){
					results.push_back(pair<int, range>(obj->id, objdis));
				}
			}
		}else{
			for(OctreeNode *c:this->children){
				c->query_within(box, results, max_dist);
			}
		}
	}
}

void OctreeNode::query_intersect(weighted_aab *box, vector<int> &results){
	if(this->intersect(*box)){
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
		root_node.update(*v);
	}
	OctreeNode *octree = new OctreeNode(root_node, 0, leaf_size);
	for(weighted_aab *v:voxels){
		octree->addObject(v);
	}
	return octree;
}

}
