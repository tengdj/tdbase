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

}

#endif /* SRC_INCLUDE_QUERY_CONTEXT_H_ */
