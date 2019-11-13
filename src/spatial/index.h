/*
 * index.h
 *
 *  Created on: Oct 30, 2019
 *      Author: teng
 */

#ifndef HISPEED_INDEX_H_
#define HISPEED_INDEX_H_

#include <stdlib.h>
#include <vector>
#include "spatial.h"

namespace hispeed{

typedef struct range{
	float low;
	float high;
} range;


class AlignedRTree{

	mbb *rtree;

	range get_range(mbb &, mbb &);


public:
	std::vector<long> get(mbb &q);

};


inline mbb get_mbb(MyMesh *mesh){
	mbb ret(mesh->bbMin, mesh->bbMax);
	return ret;
}

}
#endif /* HISPEED_INDEX_H_ */
