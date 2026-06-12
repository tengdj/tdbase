/*
 * mygpu.h
 *
 *  Created on: Dec 9, 2019
 *      Author: teng
 */

#ifndef MYGPU_H_
#define MYGPU_H_

#include <vector>
#include <mutex>
#include <thread>

using namespace std;

namespace tdbase{

class gpu_info{
public:
	int device_id;
	size_t mem_size;
	bool busy;
	mutex lock;
	char *d_data = NULL;
};
void initialize();
void init_gpu(gpu_info *gpu);
void clean_gpu(gpu_info *gpu);

vector<gpu_info *> get_gpus();
void print_gpus();

}


#endif /* MYGPU_H_ */
