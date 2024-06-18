#include "index.h"

using namespace std;

namespace tdbase{


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

inline float get_min_maxdist(vector<pair<int, range>> &results){
	float minmaxdist = DBL_MAX;
	for(int i=0;i<results.size();i++){
		minmaxdist = min(minmaxdist, results[i].second.maxdist);
	}
	return minmaxdist;
}

inline float get_max_maxdist(vector<pair<int, range>> &results){
	float maxmaxdist = 0;
	for(int i=0;i<results.size();i++){
		maxmaxdist = max(maxmaxdist, results[i].second.maxdist);
	}
	return maxmaxdist;
}

inline bool add_candidate(vector<pair<int, range>> &candidates, size_t id, range dist, int k){
	// can be inserted
	if(candidates.size()<k || dist.mindist <candidates[k-1].second.maxdist){
		int inserted_loc = 0;
		for(;inserted_loc<candidates.size();inserted_loc++){
			if(candidates[inserted_loc].first==id){
				return false;
			}
			if(candidates[inserted_loc].second.maxdist>dist.maxdist){
				candidates.insert(candidates.begin()+inserted_loc, pair<int, range>(id, dist));
				break;
			}
		}
		// insert at the end
		if(inserted_loc==candidates.size()){
			candidates.push_back(pair<int, range>(id, dist));
		}
		assert(inserted_loc<candidates.size());

		int process_loc = max(inserted_loc, k-1);

		int list_size = candidates.size();
		for(int i=process_loc+1;i<list_size;){
			// not valid anymore, remove it from the candidate list
			if(candidates[i].second.mindist >= candidates[k-1].second.maxdist){
				candidates.erase(candidates.begin()+i);
				list_size--;
				continue;
			}
			i++;
		}
		return inserted_loc<k;
	}
	return false;
}

//int ct = 0;
//inline bool update_distance_list(vector<pair<int, range>> &results, int k){
//
//	//log("%d",ct++);
//	int fartherthan = 0;
//	// check if some objects can be deleted
//	int list_size = results.size();
//	for(int i=0;i<list_size&&list_size>k;){
//		// check if this candidate still valid
//		fartherthan = 0;
//		for(int j=0;j<results.size();j++){
//			if(results[i].second>=results[j].second){
//				fartherthan++;
//			}
//		}
//		// not valid anymore, remove it from the candidate list
//		if(fartherthan>=k){
//			results.erase(results.begin()+i);
//			list_size--;
//			continue;
//		}
//		i++;
//	}
//	return true;
//
//}

void OctreeNode::query_knn(weighted_aab *box, vector<pair<int, range>> &candidates, float &min_maxdist, const int k){
	range dis = distance(*box);

	// current node possibly covers nearest objects dis.mindist<min_maxdist
	// or the candidate list is not full yet
	if(dis.mindist<min_maxdist || candidates.size()<k){
		if(isLeaf){
			for(weighted_aab *obj:objectList){
				if(obj==box){// avoid self comparing
					continue;
				}
				range objdis = obj->distance(*box);
				// min_maxdist is updated
				if(add_candidate(candidates, obj->id, objdis, k) && candidates.size()>=k){
					min_maxdist = candidates[k-1].second.maxdist;
					//log("%f",min_maxdist);
				}
			}
		}else{
			for(OctreeNode *c:this->children){
				c->query_knn(box, candidates, min_maxdist, k);
			}
		}
	}
}

void OctreeNode::query_within(weighted_aab *box, vector<pair<int, range>> &results, const float threshold){
	range dis = distance(*box);

	if(dis.mindist<=threshold){
		if(isLeaf){
			for(weighted_aab *obj:objectList){
				if(obj==box){// avoid self comparing
					continue;
				}
				range objdis = obj->distance(*box);
				if(objdis.mindist<=threshold){
					// deduplicate
					bool exist = false;
					for(pair<int, range> &c:results){
						if(c.first==obj->id){
							exist = true;
							break;
						}
					}
					if(!exist){
						results.push_back(pair<int, range>(obj->id, objdis));
					}
				}
			}
		}else{
			for(OctreeNode *c:this->children){
				c->query_within(box, results, threshold);
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
