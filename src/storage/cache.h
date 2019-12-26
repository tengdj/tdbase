/*
 * cache.h
 *
 *  Created on: Dec 26, 2019
 *      Author: teng
 */

#ifndef HISPEED_CACHE_H_
#define HISPEED_CACHE_H_

#include <pthread.h>
#include <vector>
#include "../util/util.h"

using namespace std;

namespace hispeed{


class mesh_cache{
	const size_t capacity;
	pthread_mutex_t lock;
public:
	mesh_cache(size_t c){
		capacity = c;
		pthread_mutex_init(&lock, NULL);
	}


};

}



#endif /* HISPEED_CACHE_H_ */
