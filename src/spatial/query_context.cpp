/*
 * query_context.cpp
 *
 *  Created on: Sep 24, 2022
 *      Author: teng
 */

#include <boost/program_options.hpp>
#include "query_context.h"

namespace po = boost::program_options;

namespace hispeed{

int query_context::highest_lod(){
	if(lods.size()==0){
		return 0;
	}else {
		return lods[lods.size()-1];
	}
}

query_context::query_context(){
	num_compute_thread = hispeed::get_num_threads();
	num_thread = hispeed::get_num_threads();
	pthread_mutex_init(&lk, NULL);
}

void query_context::lock(){
	pthread_mutex_lock(&lk);
}
void query_context::unlock(){
	pthread_mutex_unlock(&lk);
}

void query_context::merge(query_context ctx){
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

void query_context::report(double t){

	cout<<"total:\t"<<t<<endl;
	cout<<"index:\t"<<t*index_time/overall_time<<endl;
	cout<<"decode:\t"<<t*decode_time/overall_time<<endl;
	cout<<"packing:\t"<<t*packing_time/overall_time<<endl;
	cout<<"computation:\t"<<t*computation_time/overall_time<<endl;
	cout<<"updatelist:\t"<<t*updatelist_time/overall_time<<endl;
	cout<<"other:\t"<<t*(overall_time-index_time-decode_time-packing_time-computation_time-updatelist_time)/overall_time<<endl<<endl;
	printf("analysis\t%f\t%f\t%f\n",
			(t*index_time/overall_time)/repeated_times,
			(t*decode_time/overall_time)/repeated_times,
			(t*(computation_time+packing_time+updatelist_time)/overall_time)/repeated_times);
	cout<<"decode:\t"<<decode_time<<endl;
	cout<<"packing:\t"<<packing_time<<endl;
	printf("#objects:\t%ld\n results:\t%.2f\n", obj_count, 1.0*result_count/obj_count);
}

query_context global_ctx;

query_context parse_args(int argc, char **argv){
	query_context ctx;

	int base_lod = 0;
	int lod_gap = 50;
	int top_lod = 100;

	po::options_description desc("joiner usage");
	desc.add_options()
		("help,h", "produce help message")
		("query,q", po::value<string>(&ctx.query_type),"query type can be intersect|nn|within")
		("knn", po::value<int>(&ctx.knn), "the K value for NN query")
		("tile1", po::value<string>(&ctx.tile1_path), "path to tile 1")
		("tile2", po::value<string>(&ctx.tile2_path), "path to tile 2")
		("cn", po::value<int>(&ctx.num_compute_thread), "number of threads for geometric computation for each tile")
		("threads,n", po::value<int>(&ctx.num_thread), "number of threads for processing tiles")
		("repeat,r", po::value<int>(&ctx.repeated_times), "repeat tiles")
		("max_objects1", po::value<size_t>(&ctx.max_num_objects1), "max number of objects in tile 1")
		("max_objects2", po::value<size_t>(&ctx.max_num_objects2), "max number of objects in tile 2")
		("base_lod", po::value<int>(&base_lod), "the base lod for progressive decoding polyhedral")
		("top_lod", po::value<int>(&top_lod), "the top lod for progressive decoding polyhedral")
		("lod_gap", po::value<int>(&lod_gap), "the lod gap for progressive decoding polyhedral")
		("lod", po::value<std::vector<std::string>>()->multitoken()->
				zero_tokens()->composing(), "the lods need be processed")
		("aabb", "calculate distance with aabb")
		("gpu,g", "compute with GPU")
		("verbose,v", po::value<int>(&ctx.verbose), "verbose level")
		("counter_clock,c", "is the faces recorded clock-wise or counterclock-wise")
		("multiple_mbb,m", "using shape-aware indexing with multiple MBB")
		("max_dist", po::value<double>(&ctx.max_dist), "the maximum distance for within query")
		("hausdorf_level", po::value<int>(&ctx.hausdorf_level), "0 for no hausdor, 1 for hausdorf at the mesh level, 2 for triangle level(default)")
		;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	if (vm.count("help")) {
		cout << desc << "\n";
		exit(0);
	}
	po::notify(vm);

	if(vm.count("aabb")){
		ctx.use_aabb = true;
	}
	if(vm.count("multiple_mbb")){
		ctx.use_multimbb = true;
	}
	if(vm.count("gpu")){
		ctx.use_gpu = true;
	}
	if(vm.count("counter_clock")){
		ctx.counter_clock = true;
	}

	if(ctx.query_type!="intersect"&&ctx.query_type!="nn"&&ctx.query_type!="within"){
		cout <<"error query type: "<< ctx.query_type <<endl;
		exit(0);
	}
	if(vm.count("lod")){
		for(string l:vm["lod"].as<std::vector<std::string>>()){
			ctx.lods.push_back(atoi(l.c_str()));
		}
	}else{
		for(int l=base_lod;l<=top_lod;l+=lod_gap){
			ctx.lods.push_back(l);
		}
		if(ctx.lods[ctx.lods.size()-1]<top_lod){
			ctx.lods.push_back(top_lod);
		}
	}

	return ctx;
}


}
