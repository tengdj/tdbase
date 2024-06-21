#ifndef HISPEED_GEOMETRY_H
#define HISPEED_GEOMETRY_H

#include <math.h>
#include <stdio.h>
#include <iostream>
#include <float.h>
#include "mygpu.h"
#include "util.h"
#include "pthread.h"
using namespace std;

namespace tdbase{

// useful functions

// copy
inline
void
VcV(float Vr[3], const float V[3])
{
  Vr[0] = V[0];  Vr[1] = V[1];  Vr[2] = V[2];
}

// minus
inline
void
VmV(float Vr[3], const float V1[3], const float V2[3])
{
  Vr[0] = V1[0] - V2[0];
  Vr[1] = V1[1] - V2[1];
  Vr[2] = V1[2] - V2[2];
}

// plus
inline
void
VpV(float Vr[3], const float V1[3], const float V2[3])
{
  Vr[0] = V1[0] + V2[0];
  Vr[1] = V1[1] + V2[1];
  Vr[2] = V1[2] + V2[2];
}

// plus after product
inline
void
VpVxS(float Vr[3], const float V1[3], const float V2[3], float s)
{
  Vr[0] = V1[0] + V2[0] * s;
  Vr[1] = V1[1] + V2[1] * s;
  Vr[2] = V1[2] + V2[2] * s;
}

inline
void
VcrossV(float Vr[3], const float V1[3], const float V2[3])
{
  Vr[0] = V1[1]*V2[2] - V1[2]*V2[1];
  Vr[1] = V1[2]*V2[0] - V1[0]*V2[2];
  Vr[2] = V1[0]*V2[1] - V1[1]*V2[0];
}

// dot product
inline
float
VdotV(const float V1[3], const float V2[3])
{
  return (V1[0]*V2[0] + V1[1]*V2[1] + V1[2]*V2[2]);
}

// Euclid distance
inline
float
VdistV2(const float V1[3], const float V2[3])
{
  return ( (V1[0]-V2[0]) * (V1[0]-V2[0]) +
	   (V1[1]-V2[1]) * (V1[1]-V2[1]) +
	   (V1[2]-V2[2]) * (V1[2]-V2[2]));
}


// multiple each value in V with constant s
inline
void
VxS(float Vr[3], const float V[3], float s)
{
  Vr[0] = V[0] * s;
  Vr[1] = V[1] * s;
  Vr[2] = V[2] * s;
}

inline
void
VdS(float Vr[3], const float V[3], float s)
{
	assert(s>0);
  Vr[0] = V[0] / s;
  Vr[1] = V[1] / s;
  Vr[2] = V[2] / s;
}

class result_container{
public:
	uint32_t p1;
	uint32_t p2;
	bool intersected;
	float distance; // for normal distance
	float min_dist; // for distance range
	float max_dist;
	void print(){
		cout<<"p1:\t"<<p1<<endl;
		cout<<"p2:\t"<<p2<<endl;
		cout<<"intersected:\t"<<intersected<<endl;
		cout<<"distance:\t"<<distance<<endl;
		cout<<"min_dist:\t"<<min_dist<<endl;
		cout<<"max_dist:\t"<<max_dist<<endl;
	}
} ;

class geometry_param{
public:
	int id = 0;
	uint32_t pair_num = 0;
	uint32_t element_num = 0;
	size_t element_pair_num = 0;
	float *data = NULL;
	float *hausdorff = NULL;
	// the offset and size of the computing pairs
	uint32_t *offset_size = NULL;
	result_container *results = NULL;
	void allocate_buffer(){
		data = new float[9*element_num];
		hausdorff = new float[2*element_num];
		offset_size = new uint32_t[4*pair_num];
	}
	void clear_buffer(){
		if(data){
			delete []data;
		}
		if(hausdorff){
			delete []hausdorff;
		}
		if(offset_size){
			delete []offset_size;
		}
	}
};

inline float distance(const float *p1, const float *p2){
	float cur_dist = 0;
	for(int t=0;t<3;t++){
		cur_dist += (p1[t]-p2[t])*(p1[t]-p2[t]);
	}
	return cur_dist;
}

void compute_normal(float *Norm, const float *triangle);
bool PointInTriangleCylinder(const float *point, const float *triangle);
void project_points_to_triangle_plane(const float *point, const float *triangle, float projected_point[3]);
float PointTriangleDist(const float *point, const float *triangle);
float TriDist(const float *S, const float *T);
result_container MeshDist(const float *data1, const float *data2, size_t size1, size_t size2, const float *hausdorff1 = NULL, const float *hausdorff2 = NULL);
void MeshDist_batch_gpu(gpu_info *gpu, const float *data, const uint32_t *offset_size, const float * hausdorff, result_container *result, const uint32_t pair_num, const uint32_t element_num);

bool TriInt(const float *S, const float *T);
result_container MeshInt(const float *data1, const float *data2, size_t size1, size_t size2, const float *hausdorff1 = NULL, const float *hausdorff2 = NULL);
void TriInt_batch_gpu(gpu_info *gpu, const float *data, const uint32_t *offset_size, const float *hausdorff, result_container *result, const uint32_t batch_num, const uint32_t triangle_num);

class geometry_computer{
	pthread_mutex_t gpu_lock;
	pthread_mutex_t cpu_lock;
	int max_thread_num = tdbase::get_num_threads();
	bool cpu_busy = false;
	bool gpu_busy = false;
	bool request_cpu();
	void release_cpu();
	gpu_info *request_gpu(int min_size, bool force=false);
	void release_gpu(gpu_info *info);

	char *d_cuda = NULL;
	vector<gpu_info *> gpus;

public:
	~geometry_computer();
	geometry_computer(){
		pthread_mutex_init(&cpu_lock, NULL);
		pthread_mutex_init(&gpu_lock, NULL);
	}

	bool init_gpus();
	void get_distance_gpu(geometry_param &param);
	void get_distance_cpu(geometry_param &param);
	void get_distance(geometry_param &param);

	void get_intersect_gpu(geometry_param &param);
	void get_intersect_cpu(geometry_param &param);
	void get_intersect(geometry_param &param);
	void set_thread_num(uint32_t num){
		max_thread_num = num;
	}
};


}
#endif
