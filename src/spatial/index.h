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
#include "../geometry/aab.h"

namespace hispeed{

class AlignedRTree{
	aab *rtree;
public:
	range get_range(aab &, aab &);
	std::vector<long> get(aab &q);

};


}
#endif /* HISPEED_INDEX_H_ */
