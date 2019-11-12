/*
 * SpatialJoin.h
 *
 *  Created on: Nov 11, 2019
 *      Author: teng
 */

#ifndef SPATIALJOIN_H_
#define SPATIALJOIN_H_

#include "../voxel/voxel.h"

namespace hispeed{

// type of the workers, GPU or CPU
// each worker took a batch of jobs (X*Y) from the job queue
// and conduct the join, the result is then stored to
// the target result addresses
enum Worker_Type{
	WT_GPU,
	WT_CPU
};

// size of the buffer is 1GB
const static long VOXEL_BUFFER_SIZE = 1<<30;

class SpatialJoin{

	queue<int> set1_size;
	queue<int> set2_size;
	queue<long> set1_offset;
	queue<long> set2_offset;
	queue<float *> result_addr;
	char *buffer1;
	char *buffer2;

public:
	void SpatialJoin(){
		buffer1 = new char[VOXEL_BUFFER_SIZE];
		buffer2 = new char[VOXEL_BUFFER_SIZE];
	}
	~SpatialJoin(){
		if(buffer1){
			delete buffer1;
		}
		if(buffer2){
			delete buffer2;
		}
	}

	/*
	 * the main entry function to register the join job asynchronously
	 * others can monitor the value stored in result_addr to check whether
	 * the job is done or not
	 *
	 * */
	void register_job(int size1, int size2, const char *data1, const char *data2, float *result_addr);

	/*
	 * called by a thread doing calculation with CPUs
	 * */
	void do_job_cpu();

	/*
	 * called by a thread doing calculation with GPUs
	 * */

};



}


#endif /* SPATIALJOIN_H_ */
