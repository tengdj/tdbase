/*
 *
 * with some common gpu related operations
 * */

#include <cuda.h>
#include "cuda_util.h"
#include "mygpu.h"

using namespace std;


namespace hispeed{

vector<gpu_info *> get_gpus(){
	vector<gpu_info *> gpus;
	int num_gpus = 0;
	cudaGetDeviceCount(&num_gpus);
	for (int i = 0; i < num_gpus; i++) {
		cudaDeviceProp prop;
		cudaGetDeviceProperties(&prop, i);
		gpu_info *info = new gpu_info();
		info->busy = false;
		info->mem_size = prop.totalGlobalMem/1024/1024*4/5;
		info->device_id = i;
		// we allocate 2G mem for each gpu
		if(info->mem_size>2048){
			info->mem_size = 2048;
		}
		gpus.push_back(info);
	}
	return gpus;
}

void print_gpus(){
	int num_gpus = 0;
	cudaGetDeviceCount(&num_gpus);
	for (int i = 0; i < num_gpus; i++) {
		cudaDeviceProp prop;
		cudaGetDeviceProperties(&prop, i);
		log("Device Number: %d", i);
		log("  Device name: %s", prop.name);
		log("  Memory Clock Rate (KHz): %d", prop.memoryClockRate);
		log("  Memory Bus Width (bits): %d", prop.memoryBusWidth);
		log("  Peak Memory Bandwidth (GB/s): %f",
				2.0*prop.memoryClockRate*(prop.memoryBusWidth/8)/1.0e6);
		log("  Memory size (MB): %ld\n", prop.totalGlobalMem/1024/1024);
	}
}


void init_gpu(gpu_info *gpu){
	assert(gpu);
	if(!gpu->d_data){
		struct timeval start = get_cur_time();
		CUDA_SAFE_CALL(cudaSetDevice(gpu->device_id));
		CUDA_SAFE_CALL(cudaMalloc((void **)&gpu->d_data, gpu->mem_size*1024*1024));
		assert(gpu->d_data);
		logt("%ld MB memory size is allocated for GPU %d", start, gpu->mem_size, gpu->device_id);
	}
}

void clean_gpu(gpu_info *gpu){
	if(gpu->d_data){
		cudaSetDevice(gpu->device_id);
		CUDA_SAFE_CALL(cudaFree(gpu->d_data));
		gpu->d_data = NULL;
	}
}

void initialize(){
	cuInit(0);
}

}
