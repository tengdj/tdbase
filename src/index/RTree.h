/*
 * RTree.h
 *
 *  Created on: Oct 30, 2019
 *      Author: teng
 */

#ifndef RTREE_H_
#define RTREE_H_
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <stdlib.h>

typedef CGAL::Simple_cartesian<float> RTreeKernel;
typedef RTreeKernel::Point_3 Point;


namespace hispeed{

typedef struct mbb{
	Point min;
	Point max;
} mbb;

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

}
#endif /* RTREE_H_ */
