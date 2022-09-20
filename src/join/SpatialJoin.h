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


class query_context{
public:

	//time
	double index_time = 0;
	double decode_time = 0;
	double packing_time = 0;
	double computation_time = 0;
	double updatelist_time = 0;
	double overall_time = 0;

	//parameters
	std::string query_type = "intersect";
	int knn = 1;
	double max_dist = 1000;
	int num_thread = 0;
	int num_compute_thread = 1;
	int repeated_times = 1;
	bool use_aabb = false;
	bool use_gpu = false;
	bool use_multimbb = false;
	size_t max_num_objects1 = LONG_MAX;
	size_t max_num_objects2 = LONG_MAX;
	vector<int> lods;

	int highest_lod(){
		if(lods.size()==0){
			return 0;
		}else {
			return lods[lods.size()-1];
		}
	}

	query_context(){
		num_compute_thread = hispeed::get_num_threads();
		num_thread = hispeed::get_num_threads();
	}

	void merge(query_context ctx){
		index_time += ctx.index_time;
		decode_time += ctx.decode_time;
		packing_time += ctx.packing_time;
		computation_time += ctx.computation_time;
		updatelist_time += ctx.updatelist_time;
		overall_time += ctx.overall_time;
	}

	void report(double t){

		cout<<"total:\t"<<t<<endl;
		cout<<"index:\t"<<t*index_time/overall_time<<endl;
		cout<<"decode:\t"<<t*decode_time/overall_time<<endl;
		cout<<"packing:\t"<<t*packing_time/overall_time<<endl;
		cout<<"computation:\t"<<t*computation_time/overall_time<<endl;
		cout<<"updatelist:\t"<<t*updatelist_time/overall_time<<endl;
		cout<<"other:\t"<<t*(overall_time-decode_time-computation_time-index_time)/overall_time<<endl;
		printf("analysis\t%f\t%f\t%f\n",
				(t*index_time/overall_time)/repeated_times,
				(t*decode_time/overall_time)/repeated_times,
				(t*(computation_time+packing_time+updatelist_time)/overall_time)/repeated_times);

		cout<<"decode:\t"<<decode_time<<endl;
		cout<<"packing:\t"<<packing_time<<endl;
	}

};


class voxel_pair{
public:
	Voxel *v1;
	Voxel *v2;
	range dist;
	voxel_pair(Voxel *v1, Voxel *v2, range dist){
		this->v1 = v1;
		this->v2 = v2;
		this->dist = dist;
	};
	voxel_pair(Voxel *v1, Voxel *v2){
		this->v1 = v1;
		this->v2 = v2;
	}
};

typedef struct candidate_info_{
	HiMesh_Wrapper *mesh_wrapper;
	range distance;
	vector<voxel_pair> voxel_pairs;
}candidate_info;

typedef std::pair<HiMesh_Wrapper *, vector<candidate_info>> candidate_entry;

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

	geometry_computer *computer = NULL;

	query_context global_ctx;
	pthread_mutex_t g_lock;

public:

	SpatialJoin(geometry_computer *c, query_context &ctx){
		assert(c);
		global_ctx = ctx;
		pthread_mutex_init(&g_lock, NULL);
		computer = c;
	}
	~SpatialJoin(){}
	void report_time(double t){
		global_ctx.report(t);
	}
	/*
	 *
	 * the main entry function to conduct next round of computation
	 * each object in tile1 need to compare with all objects in tile2.
	 * to avoid unnecessary computation, we need to build index on tile2.
	 * The unit for building the index for tile2 is the ABB for all or part
	 * of the surface (mostly triangle) of a polyhedron.
	 *
	 * */
	vector<candidate_entry> mbb_knn(Tile *tile1, Tile *tile2, query_context &ctx);
	vector<candidate_entry> mbb_within(Tile *tile1, Tile *tile2, query_context &ctx);
	float *calculate_distance(vector<candidate_entry> &candidates, query_context &ctx, const int lod);
	void nearest_neighbor(Tile *tile1, Tile *tile2, query_context ctx);
	void within(Tile *tile1, Tile *tile2, query_context ctx);

	vector<candidate_entry> mbb_intersect(Tile *tile1, Tile *tile2);
	void intersect(Tile *tile1, Tile *tile2, query_context ctx);

	void join(vector<pair<Tile *, Tile *>> &tile_pairs, query_context &);

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

/* some shared functions */
void report_result(size_t id1, size_t id2);
size_t get_pair_num(vector<candidate_entry> &candidates);
size_t get_candidate_num(vector<candidate_entry> &candidates);

}


#endif /* SPATIALJOIN_H_ */
