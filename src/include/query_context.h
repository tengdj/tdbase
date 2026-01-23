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

class Configuration{
public:
	//parameters
	std::string query_type = "intersect";
	std::string tile1_path;
	std::string tile2_path;
	int knn = 1;
	float within_dist = 1000.0;
	int num_thread = 0;
	int num_compute_thread = 1;
	bool use_aabb = false;
	bool use_gpu = false;
	bool print_result = false;

	size_t max_num_objects1 = LONG_MAX;
	size_t max_num_objects2 = LONG_MAX;
	int specify_object = -1;
	vector<int> lods;
	int verbose = 0;
	bool counter_clock = false;
	bool disable_byte_encoding = false;
	Configuration(){
		num_thread = tdbase::get_num_threads();
	}

	int highest_lod(){
		if(lods.size()==0){
			return 0;
		}else {
			return lods[lods.size()-1];
		}
	}
};

class candidate_entry;
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

	// process
	uint32_t cur_lod = 0;
	vector<candidate_entry *> candidates;

	// result
	size_t obj_count = 0;
	vector<pair<int, int>> results; // contain the ID of the validate pairs
	geometry_param gp;

public:
	query_context(){
		pthread_mutex_init(&lk, NULL);
	}

	~query_context(){
		gp.clear_buffer();
		gp.clear_result();
	}

	void lock(){
		pthread_mutex_lock(&lk);
	}
	void unlock(){
		pthread_mutex_unlock(&lk);
	}

	void report_result(int id1, int id2){
		results.push_back(pair<int, int>(id1,id2));
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
		results.insert(results.end(), ctx.results.begin(), ctx.results.end());
		unlock();
	}

	void print_result(){
		sort(results.begin(),results.end());
		for(auto p:results){
			cout<<p.first<<"\t"<<p.second<<endl;
		}
	}

	void elegent_print(const char *option, double value, double overall){
		cerr<<option<<":\t"<<(value>1000?(value/1000):value)<<(value>1000?" s":" ms")<<"("<<value*100/overall<<"\%)"<<endl;
	}

	void report(){
		elegent_print("total",overall_time,overall_time);
		elegent_print("index",index_time,overall_time);
		elegent_print("decode",decode_time,overall_time);
		elegent_print("prepare",packing_time,overall_time);
		elegent_print("compute",computation_time,overall_time);
		elegent_print("evaluate",updatelist_time,overall_time);
		elegent_print("other",overall_time-index_time-decode_time-packing_time-computation_time-updatelist_time,overall_time);

		fprintf(stderr, "#objects:\t%ld\n results:%ld(\t%.3f)\n", obj_count, results.size(), 1.0*results.size()/obj_count);
	}
};

extern Configuration config;

static Configuration parse_args(int argc, char **argv){
	Configuration config;

	OptionParser op("TDBase");
	op.add<Switch>("h", "help", "produce help message");
	// for execution environment
	op.add<Value<int>>("t", "threads", "number of threads", tdbase::get_num_threads(), &config.num_thread);
	op.add<Value<int>>("", "cn", "number of threads for geometric computation for each tile", 1, &config.num_compute_thread);
	op.add<Switch>("g", "gpu", "compute with GPU", &config.use_gpu);
	op.add<Value<int>>("v", "verbose", "verbose level", 0, &config.verbose);
	op.add<Switch>("p", "print_result", "print result to standard out", &config.print_result);
	// for system configuration
	op.add<Switch>("", "aabb", "calculate distance with aabb", &config.use_aabb);
	op.add<Switch>("c", "counter_clock", "is the faces recorded clock-wise or counterclock-wise", &config.counter_clock);
	op.add<Switch>("", "disable_byte_encoding", "using the raw hausdorff instead of the byte encoded ones", &config.disable_byte_encoding);
	auto lod_options = op.add<Value<int>>("l", "lods", "the lods that needs be processed");

	// for input data
	op.add<Value<string>, Attribute::required>("", "tile1", "path to tile 1", "",  & config.tile1_path);
	op.add<Value<string>>("", "tile2", "path to tile 2", "",  & config.tile2_path);
	op.add<Value<size_t>>("", "max_objects1", "max number of objects in tile 1", LONG_MAX, &config.max_num_objects1);
	op.add<Value<size_t>>("", "max_objects2", "max number of objects in tile 2", LONG_MAX, &config.max_num_objects2);
	op.add<Value<int>>("", "specify_object", "specify a single object in tile 1 for processing", -1, &config.specify_object);
	// for query 
	op.add<Value<string>, Attribute::required>("q", "query", "query types: intersect|nn|within", "", & config.query_type);
	op.add<Value<int>>("k", "knn", "the K value for NN query", 1, &config.knn);
	op.add<Value<float>>("w", "within_dist", "the maximum distance for within query", 1000.0, &config.within_dist);

	op.parse(argc, argv);

// post process

	if(config.query_type!="intersect"&&config.query_type!="nn"&&config.query_type!="within"){
		cout <<"error query type: "<< config.query_type <<endl;
		exit(0);
	}

	if(lod_options->count()>0){
		for (int i = 0; i < lod_options->count(); i++) {
			config.lods.push_back(lod_options->value(i));
		}
	}else{
		for(int l=20;l<=100;l+=20){
			config.lods.push_back(l);
		}
	}
	sort(config.lods.begin(), config.lods.end());
	unique(config.lods.begin(), config.lods.end());

	return config;
}

}

#endif /* SRC_INCLUDE_QUERY_CONTEXT_H_ */
