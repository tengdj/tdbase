/*
 * query_context.h
 *
 *  Created on: Sep 21, 2022
 *      Author: teng
 */

#ifndef SRC_INCLUDE_QUERY_CONTEXT_H_
#define SRC_INCLUDE_QUERY_CONTEXT_H_

#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <iostream>
#include "util.h"
#include "geometry.h"

using namespace std;

namespace hispeed{

class Tile;
class query_context{
public:
	pthread_mutex_t lk;

	//time
	double index_time = 0;
	double decode_time = 0;
	double packing_time = 0;
	double computation_time = 0;
	double updatelist_time = 0;
	double overall_time = 0;

	//parameters
	std::string query_type = "intersect";
	std::string tile1_path;
	std::string tile2_path;
	int knn = 1;
	double max_dist = 1000;
	int num_thread = 0;
	int num_compute_thread = 1;
	int repeated_times = 1;
	bool use_aabb = false;
	bool use_gpu = false;
	bool use_multimbb = false;
	int hausdorf_level = 2; // 0 for no hausdor, 1 for hausdorf at the mesh level, 2 for triangle level
	size_t max_num_objects1 = LONG_MAX;
	size_t max_num_objects2 = LONG_MAX;
	vector<int> lods;
	int verbose = 0;
	bool counter_clock = false;

	uint cur_lod = 0;
	Tile *tile1 = NULL;
	Tile *tile2 = NULL;

	// result
	size_t obj_count = 0;
	size_t result_count = 0;
	result_container *results = NULL;

	int highest_lod();
	query_context();
	void merge(query_context ctx);
	void report(double t);
	void lock();
	void unlock();
};

extern query_context global_ctx;
extern query_context parse_args(int argc, char **argv);

}

#endif /* SRC_INCLUDE_QUERY_CONTEXT_H_ */
