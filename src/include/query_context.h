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
#include "popl.h"

using namespace std;
using namespace popl;

namespace tdbase{

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
	float within_dist = 1000.0;
	int num_thread = 0;
	int num_compute_thread = 1;
	int repeated_times = 1;
	bool use_aabb = false;
	bool use_gpu = false;
	bool print_result = false;
	int hausdorf_level = 2; // 0 for no hausdorff, 1 for hausdorff at the mesh level, 2 for triangle level hausdorff
	size_t max_num_objects1 = LONG_MAX;
	size_t max_num_objects2 = LONG_MAX;
	vector<int> lods;
	int verbose = 0;
	bool counter_clock = false;
	bool disable_byte_encoding = false;

	uint32_t cur_lod = 0;
	Tile *tile1 = NULL;
	Tile *tile2 = NULL;

	// result
	size_t obj_count = 0;
	size_t result_count = 0;
	result_container *results = NULL;

	query_context(){
		num_thread = tdbase::get_num_threads();
		pthread_mutex_init(&lk, NULL);
	}

	void lock(){
		pthread_mutex_lock(&lk);
	}
	void unlock(){
		pthread_mutex_unlock(&lk);
	}

	int highest_lod(){
		if(lods.size()==0){
			return 0;
		}else {
			return lods[lods.size()-1];
		}
	}

	void merge(query_context &ctx){
		lock();
		index_time += ctx.index_time;
		decode_time += ctx.decode_time;
		packing_time += ctx.packing_time;
		computation_time += ctx.computation_time;
		updatelist_time += ctx.updatelist_time;
		overall_time += ctx.overall_time;
		obj_count += ctx.obj_count;
		result_count += ctx.result_count;
		unlock();
	}

	void report(double t){

		cerr<<"total:\t"<<t<<endl;
		cerr<<"index:\t"<<t*index_time/overall_time<<endl;
		cerr<<"decode:\t"<<t*decode_time/overall_time<<endl;
		cerr<<"packing:\t"<<t*packing_time/overall_time<<endl;
		cerr<<"computation:\t"<<t*computation_time/overall_time<<endl;
		cerr<<"updatelist:\t"<<t*updatelist_time/overall_time<<endl;
		cerr<<"other:\t"<<t*(overall_time-index_time-decode_time-packing_time-computation_time-updatelist_time)/overall_time<<endl<<endl;
		fprintf(stderr, "analysis\t%f\t%f\t%f\n",
				(t*index_time/overall_time)/repeated_times,
				(t*(decode_time+packing_time+updatelist_time)/overall_time)/repeated_times,
				(t*computation_time/overall_time)/repeated_times);
		cerr<<"decode:\t"<<decode_time<<endl;
		cerr<<"packing:\t"<<packing_time<<endl;
		fprintf(stderr, "#objects:\t%ld\n results:%ld(\t%.3f)\n", obj_count, result_count, 1.0*result_count/obj_count);

		fprintf(stderr, "%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
				t*index_time/overall_time,
				t*decode_time/overall_time,
				t*packing_time/overall_time,
				t*computation_time/overall_time,
				t*updatelist_time/overall_time,
				t*(overall_time-index_time-decode_time-packing_time-computation_time-updatelist_time)/overall_time,
				t);
	}
};

extern query_context global_ctx;

static query_context parse_args(int argc, char **argv){
	query_context ctx;

	OptionParser op("TDBase");
	op.add<Switch>("h", "help", "produce help message");
	// for execution environment
	op.add<Value<int>>("t", "threads", "number of threads", tdbase::get_num_threads(), &ctx.num_thread);
	op.add<Value<int>>("", "cn", "number of threads for geometric computation for each tile", 1, &ctx.num_compute_thread);
	op.add<Switch>("g", "gpu", "compute with GPU", &ctx.use_gpu);
	op.add<Value<int>>("v", "verbose", "verbose level", 0, &ctx.verbose);
	op.add<Switch>("", "print_result", "print result to standard out", &ctx.print_result);
	// for system configuration
	op.add<Switch>("", "aabb", "calculate distance with aabb", &ctx.use_aabb);
	op.add<Switch>("c", "counter_clock", "is the faces recorded clock-wise or counterclock-wise", &ctx.counter_clock);
	op.add<Switch>("", "disable_byte_encoding", "using the raw hausdorff instead of the byte encoded ones", &ctx.disable_byte_encoding);
	op.add<Value<int>>("", "hausdorff_level", "0 for no hausdorff, 1 for hausdorff at the mesh level, 2 for triangle level", 2, &ctx.hausdorf_level);
	auto lod_options = op.add<Value<int>>("l", "lods", "the lods that needs be processed");

	// for input data
	op.add<Value<string>, Attribute::required>("", "tile1", "path to tile 1", "",  & ctx.tile1_path);
	op.add<Value<string>>("", "tile2", "path to tile 2", "",  & ctx.tile2_path);
	op.add<Value<size_t>>("", "max_objects1", "max number of objects in tile 1", LONG_MAX, &ctx.max_num_objects1);
	op.add<Value<size_t>>("", "max_objects2", "max number of objects in tile 2", LONG_MAX, &ctx.max_num_objects2);
	// for query 
	op.add<Value<string>, Attribute::required>("q", "query", "query types: intersect|nn|within", "", & ctx.query_type);
	op.add<Value<int>>("k", "knn", "the K value for NN query", 1, &ctx.knn);
	op.add<Value<float>>("w", "within_dist", "the maximum distance for within query", 1000.0, &ctx.within_dist);

	op.parse(argc, argv);

// post process
	assert(ctx.hausdorf_level>=0 && ctx.hausdorf_level<=2);

	if(ctx.query_type!="intersect"&&ctx.query_type!="nn"&&ctx.query_type!="within"){
		cout <<"error query type: "<< ctx.query_type <<endl;
		exit(0);
	}
	if(lod_options->count()>0){
		for (int i = 0; i < lod_options->count(); i++) {
			ctx.lods.push_back(lod_options->value(i));
		}
	}else{
		for(int l=20;l<=100;l+=20){
			ctx.lods.push_back(l);
		}
	}
	sort(ctx.lods.begin(), ctx.lods.end());
	unique(ctx.lods.begin(), ctx.lods.end());

	global_ctx.num_thread = min(global_ctx.num_thread, global_ctx.repeated_times);
	return ctx;
}

}

#endif /* SRC_INCLUDE_QUERY_CONTEXT_H_ */
