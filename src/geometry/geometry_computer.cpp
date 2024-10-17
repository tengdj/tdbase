/*
 * geometry_computer.cpp
 *
 *  Created on: Dec 9, 2019
 *      Author: teng
 */

#include "geometry.h"
#include "mygpu.h"
#include "query_context.h"


namespace tdbase{


void *MeshDist_unit(void *params_void){
	geometry_param *param = (geometry_param *)params_void;
	for(int i=0;i<param->pair_num;i++){
		param->results[i] = MeshDist(param->data+param->offset_size[4*i]*9,
									   param->data+param->offset_size[4*i+2]*9,
									   param->offset_size[4*i+1],
									   param->offset_size[4*i+3],
									   param->hausdorff+param->offset_size[4*i]*2,
									   param->hausdorff+param->offset_size[4*i+2]*2);

	}
	return NULL;
}

bool geometry_computer::request_cpu(){
	if(cpu_busy==false){
		pthread_mutex_lock(&cpu_lock);
		if(cpu_busy==false){
			cpu_busy = true;
			return true;
		}else{
			pthread_mutex_unlock(&cpu_lock);
		}
	}
	return false;
}
void geometry_computer::release_cpu(){
	assert(cpu_busy);
	cpu_busy = false;
	pthread_mutex_unlock(&cpu_lock);
}

#ifdef USE_GPU

gpu_info *geometry_computer::request_gpu(int min_size, bool force){
	assert(gpus.size()>=1);
	do{
		for(gpu_info *info:gpus){
			if(!info->busy&&info->mem_size>min_size+1){
				pthread_mutex_lock(&info->lock);
				if(!info->busy){
					info->busy = true;
					// will be released in release_gpu() function
					return info;
				}
				pthread_mutex_unlock(&info->lock);
			}
		}
	}while(force);
	return NULL;
}

void geometry_computer::release_gpu(gpu_info *info){
	assert(info->busy);
	pthread_mutex_unlock(&info->lock);
	info->busy = false;
}

bool geometry_computer::init_gpus(){
	gpus = get_gpus();
	for(gpu_info *info:gpus){
		init_gpu(info);
	}
	return true;
}

void geometry_computer::get_distance_gpu(geometry_param &cc){
	gpu_info *gpu = request_gpu(cc.element_num*11*sizeof(float)/1024/1024, true);
	assert(gpu);
	log("GPU %d started to get distance", gpu->device_id);
	tdbase::MeshDist_batch_gpu(gpu, cc.data, cc.offset_size, cc.hausdorff, cc.results, cc.pair_num, cc.element_num);
	release_gpu(gpu);
}


void geometry_computer::get_intersect_gpu(geometry_param &cc){
	gpu_info *gpu = request_gpu(cc.element_num*9*sizeof(float)/1024/1024, true);
	assert(gpu);
	log("GPU %d started to check intersect", gpu->device_id);
	tdbase::TriInt_batch_gpu(gpu, cc.data, cc.offset_size, cc.hausdorff, cc.results, cc.pair_num, cc.element_num);
	release_gpu(gpu);
}

#endif

geometry_computer::~geometry_computer(){
#ifdef USE_GPU
	for(gpu_info *info:gpus){
		clean_gpu(info);
		delete info;
	}
#endif
}
void geometry_computer::get_distance_cpu(geometry_param &cc){

	int each_thread = cc.pair_num/max_thread_num;
	int thread_num = max_thread_num;
	if(each_thread==0){
		thread_num = cc.pair_num;
	}
	each_thread++;
	pthread_t threads[thread_num];
	geometry_param params[thread_num];

	int i=0;
	for(;i<thread_num;i++){
		int start = each_thread*i;
		if(start>=cc.pair_num){
			break;
		}

		params[i] = cc;
		params[i].pair_num = min(each_thread, (int)cc.pair_num-start);
		params[i].offset_size = cc.offset_size+start*4;
		params[i].data = cc.data;
		params[i].hausdorff = cc.hausdorff;
		params[i].id = i+1;
		params[i].results = cc.results+start;
		pthread_create(&threads[i], NULL, MeshDist_unit, (void *)&params[i]);
	}
	if(max_thread_num>1){
		log("%d threads started to get distance", thread_num);
	}
	for(; i > 0; i--){
		void *status;
		pthread_join(threads[i-1], &status);
	}
}



void geometry_computer::get_distance(geometry_param &cc){
	if(gpus.size()>0){
#ifdef USE_GPU
		get_distance_gpu(cc);
#endif
	}else{
		get_distance_cpu(cc);
	}
//	for(int i=0;i<cc.pair_num;i++){
//		cc.results[i].print();
//	}
}

void *TriInt_unit(void *params_void){
	geometry_param *param = (geometry_param *)params_void;
	for(int i=0;i<param->pair_num;i++){
		param->results[i] = MeshInt(param->data+param->offset_size[4*i]*9,
									    param->data+param->offset_size[4*i+2]*9,
									    param->offset_size[4*i+1],
									    param->offset_size[4*i+3],
										param->hausdorff+param->offset_size[4*i]*2,
										param->hausdorff+param->offset_size[4*i+2]*2);
	}
	return NULL;
}

void geometry_computer::get_intersect_cpu(geometry_param &cc){

	pthread_t threads[max_thread_num];
	geometry_param params[max_thread_num];
	int each_thread = cc.pair_num/max_thread_num+1;
	int tnum = 0;
	for(tnum=0;tnum<max_thread_num;tnum++){
		int start = each_thread*tnum;
		if(start>=cc.pair_num){
			break;
		}
		params[tnum] = cc;
		params[tnum].id = tnum+1;
		params[tnum].pair_num = min(each_thread, (int)cc.pair_num-start);
		params[tnum].offset_size = cc.offset_size+start*4;
		params[tnum].data = cc.data;
		params[tnum].hausdorff = cc.hausdorff;
		params[tnum].results = cc.results+start;
		pthread_create(&threads[tnum], NULL, TriInt_unit, (void *)&params[tnum]);
	}
	log("%d threads started", max_thread_num);
	for(int i = 0; i < tnum; i++){
		void *status;
		pthread_join(threads[i], &status);
	}
}

void geometry_computer::get_intersect(geometry_param &cc){

	if(gpus.size()>0){
#ifdef USE_GPU
		get_intersect_gpu(cc);
#endif
	}else{
		get_intersect_cpu(cc);
	}
}

}

