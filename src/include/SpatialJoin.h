/*
 * SpatialJoin.h
 *
 *  Created on: Nov 11, 2019
 *      Author: teng
 */

#ifndef SPATIALJOIN_H_
#define SPATIALJOIN_H_

#include <queue>
#include "query_context.h"
#include "aab.h"
#include "tile.h"
#include "geometry.h"
#include "candidate.h"
#include "himesh.h"
using namespace std;

namespace tdbase{

// type of the workers, GPU or CPU
// each worker took a batch of jobs (X*Y) from the job queue
// and conduct the join, the result is then stored to
// the target result addresses
enum Worker_Type{
	WT_GPU,
	WT_CPU
};

enum Join_Type{
	JT_intersect,
	JT_distance,
	JT_nearest
};

// size of the buffer is 1GB
const static long VOXEL_BUFFER_SIZE = 1<<30;

/* todo: need be updated
 * all computations will be aligned into computation units
 * of N*N. For instance, after checking the index, the
 * nearest neighbor of object a is b or c is not decided.
 * We further decode the polyhedron a, b and c if they
 * are not decoded yet. The edges and surfaces are
 * decoded and cached. Then the computation across those
 * segments and triangles are organized as many N*N
 * computing units as possible. Notice that padding may be needed to
 * align the computing units. then space in buffer is claimed
 * to store the data of those computing units. the true computation
 * is done by the GPU or CPU, and results will be copied into the result_addr
 * Corresponding to the buffer space claimed.
*/

class SpatialJoin{
protected:
	geometry_computer *computer = NULL;
	// five steps of query
	virtual void index_retrieval(Tile *tile1, Tile *tile2, query_context &ctx) = 0;
	void decode_data(query_context &ctx);
	void packing_data(query_context &ctx);
	virtual void geometric_computation(query_context &ctx)=0;
	virtual void evaluate_candidate_lists(query_context &ctx)=0;
public:
	SpatialJoin();
	virtual ~SpatialJoin();
	void join(Tile *, Tile *);
};

class DistanceJoin:public SpatialJoin{
protected:
	range update_voxel_pair_list(vector<voxel_pair> &voxel_pairs, double minmaxdist, bool keep_empty=true);
	void update_distance_ranges(query_context &ctx);
	void geometric_computation( query_context &ctx);
public:
	DistanceJoin():SpatialJoin(){}
};

class KNNJoin:public DistanceJoin{
protected:
	void index_retrieval(Tile *tile1, Tile *tile2, query_context &ctx);
	void evaluate_candidate_lists(query_context &ctx);
public:
	KNNJoin():DistanceJoin(){}
};

class DWithinJoin:public DistanceJoin{
protected:
	void index_retrieval(Tile *tile1, Tile *tile2, query_context &ctx);
	void evaluate_candidate_lists(query_context &ctx);
public:
	DWithinJoin():DistanceJoin(){}
};

class IntersectJoin:public SpatialJoin{
protected:
	void index_retrieval(Tile *tile1, Tile *tile2, query_context &ctx);
	void geometric_computation(query_context &ctx);
	void evaluate_candidate_lists(query_context &ctx);
public:
	IntersectJoin():SpatialJoin(){}
};

}


#endif /* SPATIALJOIN_H_ */
