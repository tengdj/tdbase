/*
 * SpatialJoin.h
 *
 *  Created on: Nov 11, 2019
 *      Author: teng
 */

#ifndef SPATIALJOIN_H_
#define SPATIALJOIN_H_

#include "../storage/tile.h"
#include "../geometry/geometry.h"
#include <queue>

using namespace std;

namespace hispeed{

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

class SpatialJoin{

	// tiles with data for join
	Tile *tile1;
	Tile *tile2;

	geometry_computer *computer = NULL;

	enum Join_Type type = JT_nearest;
	// cursor of the object in tile 1
	// waiting for processing
	long cursor = 0;
	/*
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
	float *result_addr = NULL;
	float *buffer = NULL;

	// sign for completeness
	bool complete = false;

	float *d_data = NULL;
public:

	SpatialJoin(Tile *t1, Tile *t2, geometry_computer *c){
		assert(t1!=NULL&&t2!=NULL&&c);
		tile1 = t1;
		tile2 = t2;
		computer = c;
	}
	~SpatialJoin(){
		if(buffer){
			delete buffer;
		}
	}

	inline bool is_self_join(){
		return tile1==tile2;
	};

	/*
	 *
	 * the main entry function to conduct next round of computation
	 * each object in tile1 need to compare with all objects in tile2.
	 * to avoid unnecessary computation, we need to build index on tile2.
	 * The unit for building the index for tile2 is the ABB for all or part
	 * of the surface (mostly triangle) of a polyhedron.
	 *
	 * */
	void nearest_neighbor(bool with_gpu);
	void intersect(bool with_gpu);


	/*
	 *
	 * go check the index
	 *
	 * */
	void check_index();

	// register job to gpu
	// worker can register work to gpu
	// if the GPU is idle and queue is not full
	// otherwise do it locally with CPU
	float *register_computation(char *data, int num_cu);

	/*
	 * do the geometry computation in a batch with GPU
	 *
	 * */
	void compute_gpu();

	/*
	 * do the geometry computation in a batch with CPU
	 *
	 * */
	void compute_cpu();


};



}


#endif /* SPATIALJOIN_H_ */
