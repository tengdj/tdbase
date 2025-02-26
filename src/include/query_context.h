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
#include <boost/program_options.hpp>

#include "util.h"
#include "geometry.h"

using namespace std;

namespace po = boost::program_options;

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
	double within_dist = 1000;
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

	po::options_description desc("joiner usage");
	desc.add_options()
		("help,h", "produce help message")
		("aabb", "calculate distance with aabb")
		("gpu,g", "compute with GPU")
		("counter_clock,c", "is the faces recorded clock-wise or counterclock-wise")
		("disable_byte_encoding", "using the raw hausdorff instead of the byte encoded ones")

		// for data
		("tile1", po::value<string>(&ctx.tile1_path), "path to tile 1")
		("tile2", po::value<string>(&ctx.tile2_path), "path to tile 2")
		("repeat,r", po::value<int>(&ctx.repeated_times), "repeat tiles")
		("max_objects1", po::value<size_t>(&ctx.max_num_objects1), "max number of objects in tile 1")
		("max_objects2", po::value<size_t>(&ctx.max_num_objects2), "max number of objects in tile 2")

		// query setup
		("query,q", po::value<string>(&ctx.query_type),"query type can be intersect|nn|within")
		("knn", po::value<int>(&ctx.knn), "the K value for NN query")
		("within_dist", po::value<double>(&ctx.within_dist), "the maximum distance for within query")
		("lod", po::value<std::vector<std::string>>()->multitoken()->zero_tokens()->composing(), "the lods need be processed")
		("hausdorf_level", po::value<int>(&ctx.hausdorf_level), "0 for no hausdorff, 1 for hausdorff at the mesh level, 2 for triangle level(default)")

		// execution setup
		("cn", po::value<int>(&ctx.num_compute_thread), "number of threads for geometric computation for each tile")
		("threads,n", po::value<int>(&ctx.num_thread), "number of threads for processing tiles")
		("verbose,v", po::value<int>(&ctx.verbose), "verbose level")		
		("print_result", "print result to standard out")
		;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	if (vm.count("help")) {
		cout << desc << "\n";
		exit(0);
	}
	po::notify(vm);

	ctx.use_aabb = vm.count("aabb");
	ctx.use_gpu = vm.count("gpu");
	ctx.counter_clock = vm.count("counter_clock");
	ctx.disable_byte_encoding = vm.count("disable_byte_encoding");
	ctx.print_result = vm.count("print_result");

	assert(ctx.hausdorf_level>=0 && ctx.hausdorf_level<=2);

	if(ctx.query_type!="intersect"&&ctx.query_type!="nn"&&ctx.query_type!="within"){
		cout <<"error query type: "<< ctx.query_type <<endl;
		exit(0);
	}
	if(vm.count("lod")){
		for(string l:vm["lod"].as<std::vector<std::string>>()){
			ctx.lods.push_back(atoi(l.c_str()));
		}
	}else{
		for(int l=20;l<=100;l+=20){
			ctx.lods.push_back(l);
		}
	}

	global_ctx.num_thread = min(global_ctx.num_thread, global_ctx.repeated_times);

	return ctx;
}

}

#endif /* SRC_INCLUDE_QUERY_CONTEXT_H_ */
