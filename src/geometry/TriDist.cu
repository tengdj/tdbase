/*************************************************************************\

  Copyright 1999 The University of North Carolina at Chapel Hill.
  All Rights Reserved.

  Permission to use, copy, modify and distribute this software and its
  documentation for educational, research and non-profit purposes, without
  fee, and without a written agreement is hereby granted, provided that the
  above copyright notice and the following three paragraphs appear in all
  copies.

  IN NO EVENT SHALL THE UNIVERSITY OF NORTH CAROLINA AT CHAPEL HILL BE
  LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR
  CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE
  USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY
  OF NORTH CAROLINA HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH
  DAMAGES.

  THE UNIVERSITY OF NORTH CAROLINA SPECIFICALLY DISCLAIM ANY
  WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE SOFTWARE
  PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
  NORTH CAROLINA HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT,
  UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

  The authors may be contacted via:

  US Mail:             E. Larsen
                       Department of Computer Science
                       Sitterson Hall, CB #3175
                       University of N. Carolina
                       Chapel Hill, NC 27599-3175

  Phone:               (919)962-1749

  EMail:               geom@cs.unc.edu


\**************************************************************************/

//--------------------------------------------------------------------------
// File:   TriDist.cpp
// Author: Eric Larsen
// Description:
// contains SegPoints() for finding closest points on a pair of line
// segments and TriDist() for finding closest points on a pair of triangles
//--------------------------------------------------------------------------

#include "TriDist.h"
#include <pthread.h>

#define CUDA_ERROR_CHECK

#define CudaSafeCall( err ) __cudaSafeCall( err, __FILE__, __LINE__ )
#define CudaCheckError()    __cudaCheckError( __FILE__, __LINE__ )

namespace hispeed{

inline void __cudaSafeCall( cudaError err)
{
#ifdef CUDA_ERROR_CHECK
    if ( cudaSuccess != err )
    {
        fprintf( stderr, "cudaSafeCall() failed: %s\n",
                 cudaGetErrorString( err ) );
        exit( -1 );
    }
#endif

    return;
}

// copy
__device__
inline void
VcV(float Vr[3], const float V[3])
{
  Vr[0] = V[0];  Vr[1] = V[1];  Vr[2] = V[2];
}

// minus
__device__
inline void
VmV(float Vr[3], const float V1[3], const float V2[3])
{
  Vr[0] = V1[0] - V2[0];
  Vr[1] = V1[1] - V2[1];
  Vr[2] = V1[2] - V2[2];
}

// plus
__device__
inline void
VpV(float Vr[3], const float V1[3], const float V2[3])
{
  Vr[0] = V1[0] + V2[0];
  Vr[1] = V1[1] + V2[1];
  Vr[2] = V1[2] + V2[2];
}

// plus after product
__device__
inline void
VpVxS(float Vr[3], const float V1[3], const float V2[3], float s)
{
  Vr[0] = V1[0] + V2[0] * s;
  Vr[1] = V1[1] + V2[1] * s;
  Vr[2] = V1[2] + V2[2] * s;
}

// dot product
__device__
inline float
VdotV(const float V1[3], const float V2[3])
{
  return (V1[0]*V2[0] + V1[1]*V2[1] + V1[2]*V2[2]);
}

__device__
inline float SegDist_kernel(const float *S, const float *T,
							const float *A, const float *B)
{

	float t = 0.0, u = 0.0, dist = 0.0, t1 = 0.0;
	float ST[3]; // temporary vector S->T
	VmV(ST,T,S);
	float A_dot_A = VdotV(A,A);
	float B_dot_B = VdotV(B,B);
	if(A_dot_A==0||B_dot_B==0){
		return DBL_MAX;
	}
	float A_dot_B = VdotV(A,B);
	float A_dot_ST = VdotV(A,ST);
	float B_dot_ST = VdotV(B,ST);

	// t parameterizes ray P,A
	// u parameterizes ray Q,B

	// compute t for the closest point on ray P,A to
	// ray Q,B

	float denom = A_dot_A*B_dot_B - A_dot_B*A_dot_B;
	if(denom == 0){
		t = 0;
	}else{
		t = (A_dot_ST*B_dot_B - B_dot_ST*A_dot_B) / denom;
	}

	// find u for point on ray Q,B closest to point at t
	// B_dot_B can never be 0
	u = (t*A_dot_B - B_dot_ST)/B_dot_B;
	// if u is on segment Q,B, t and u correspond to
	// closest points, otherwise, recompute and
	// clamp t

	if (u <= 0) {
		u = 0;
	} else if (u >= 1) {
		u = 1;
	}

	t = (A_dot_B*u+A_dot_ST)/A_dot_A;

	if(t<=0){
		t = 0;
	} else if(t >= 1){
		t = 1;
	}

	t1 = A[0]*t-ST[0]-B[0]*u;
	dist += t1*t1;
	t1 = A[1]*t-ST[1]-B[1]*u;
	dist += t1*t1;
	t1 = A[2]*t-ST[2]-B[2]*u;
	dist += t1*t1;

	return dist;
}


//	float idf = (float)id;
//	float AA[3]={1,1,1},BB[3]={1,1,1},
//			SS[6]={2+idf,2+idf,2+idf,3+idf,3+idf,3+idf},
//			TT[6]={0+idf,0+idf,0+idf,1+idf,1+idf,1+idf};
//	float dd = SegDist_kernel(SS, TT, AA, BB);


__global__
void SegDist_cuda(const float *S, const float *T,
				  const float *A, const float *B,
				  float *dist){
	int x = blockIdx.x;
	int y = threadIdx.x;
	int id = x*blockDim.x+y;
	const float *cur_S = S+x*6;
	const float *cur_T = T+y*6;
	const float *cur_A = A+x*3;
	const float *cur_B = B+y*3;
	float dd = SegDist_kernel(cur_S, cur_T, cur_A, cur_B);

	if(dist[id]>dd){
		dist[id] = dd;
	}
}


__global__
void SegDist_cuda_new(const float *S, const float *T,
				  	  const float *A, const float *B,
				  	  float *dist, int voxel_size){
	// computing which segment in set1 of voxel with
	// voxel_id
	int voxel_id = blockIdx.x;
	int segment_id = threadIdx.x;
	// update the pointers for current thread
	const float *cur_S = S+6*voxel_id*400+6*segment_id;
	const float *cur_T = T+6*voxel_id*400;
	const float *cur_A = A+3*voxel_id*400+3*segment_id;
	const float *cur_B = B+3*voxel_id*400;

	float min_dist = DBL_MAX;
	// go over all the segments in set2
	for(int i = 0;i<voxel_size;i++){
		float dd = SegDist_kernel(cur_S, cur_T, cur_A, cur_B);
		if(min_dist>dd){
			min_dist = dd;
		}
		cur_T += 6;
		cur_B += 3;
	}
	// initialize the minimum distance
	if(segment_id == 0){
		dist[voxel_id] = min_dist;
	}
	// update the minimum distance
	for(int i=1;i<voxel_size;i++){
		if(i==segment_id){
			if(dist[voxel_id]>min_dist){
				dist[voxel_id] = min_dist;
			}
		}
	}
}

__global__
void get_max(float *d, float *max_d, int batch)
{
	int id = threadIdx.x;
	float min_val = DBL_MAX;
	float *cur_d = d+batch*id;
	for(int i=0;i<batch;i++){
		if(cur_d[i]<min_val){
			min_val = cur_d[i];
		}
	}
	max_d[id] = min_val;
}

__global__
void get_vector_kernel(float *S, float *A){
	int id = blockIdx.x*blockDim.x+threadIdx.x;
	float *cur_A = A+id*3;
	VmV(cur_A, S+id*6+3, S+id*6);
}


// the computing capacity of the gpu
int max_len_x = 200;
int max_len_y = 512;

float SegDist_batch_gpu(const float *S, const float *T, int size1, int size2){

	struct timeval start = get_cur_time();
	float *d_S, *d_T, *d_A, *d_B;
	float *d_dist;

	int len_x = min(max_len_x, size1);
	int len_y = min(max_len_y, size2);
	float min_dist = DBL_MAX;

	// copy data into device
	cudaMalloc(&d_S, size1*6*sizeof(float));
	cudaMalloc(&d_T, size2*6*sizeof(float));
	cudaMemcpy(d_S, S, size1*6*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_T, T, size2*6*sizeof(float), cudaMemcpyHostToDevice);

	// some temporary space for computation
	cudaMalloc(&d_A, size1*3*sizeof(float));
	cudaMalloc(&d_B, size2*3*sizeof(float));

	// compute the vectors of segments in S and T, save to A and B
	for(int i=0;i<size1;){
		float *cur_S = d_S+i*6;
		float *cur_A = d_A+i*3;

		int b_num;
		int t_num;
		if(size1>i+max_len_x*max_len_y){
			b_num = max_len_x;
			t_num = max_len_y;
			i += max_len_x*max_len_y;
		}else if(size1 < i+max_len_y){
			b_num = 1;
			t_num = size1-i;
			i = size1;
		}else{
			b_num = (size1-i)/max_len_y;
			t_num = max_len_y;
			i += b_num*max_len_y;
		}
		get_vector_kernel<<<b_num, t_num>>>(cur_S, cur_A);
	}
	for(int i=0;i<size2;){
		float *cur_T = d_T+i*6;
		float *cur_B = d_B+i*3;

		int b_num;
		int t_num;
		if(size2>i+max_len_x*max_len_y){
			b_num = max_len_x;
			t_num = max_len_y;
			i += max_len_x*max_len_y;
		}else if(size2 < i+max_len_y){
			b_num = 1;
			t_num = size2-i;
			i = size2;
		}else{
			b_num = (size2-i)/max_len_y;
			t_num = max_len_y;
			i += b_num*max_len_y;
		}
		get_vector_kernel<<<b_num, t_num>>>(cur_T, cur_B);
	}
	cudaDeviceSynchronize();

	// space for the distances got from gpu
	cudaMalloc(&d_dist, len_x*len_y*sizeof(float));
	float *dist = new float[len_x*len_y];
	for(int i=0;i<len_x*len_y;i++){
		dist[i] = DBL_MAX;
	}
	cudaMemcpy(d_dist, dist, len_x*len_y*sizeof(float), cudaMemcpyHostToDevice);
	cout<<"preprocessing data takes "<<get_time_elapsed(start)<<" ms"<<endl;

	start = get_cur_time();
	const float *cur_S, *cur_T, *cur_A, *cur_B;
	int times = 0;
	for(int i=0;i<size1;i+=len_x){
		for(int j=0;j<size2;j+=len_y){
			times++;
			cur_S = d_S+i*6;
			cur_T = d_T+j*6;
			cur_A = d_A+i*3;
			cur_B = d_B+j*3;
			int tsize1 = min(len_x, size1-i);
			int tsize2 = min(len_y, size2-j);
			SegDist_cuda<<<tsize1, tsize2>>>(cur_S, cur_T, cur_A, cur_B, d_dist);
		}
	}
	cudaDeviceSynchronize();
	double time_elapsed = get_time_elapsed(start);

	cout<<"run "<<times<<" rounds in "<<time_elapsed<<" ms, each round takes "<<time_elapsed/times<<" ms "<<endl;

	start = get_cur_time();
	cudaMemcpy(dist, d_dist, len_x*len_y*sizeof(float), cudaMemcpyDeviceToHost);
	for(int i=0;i<len_x*len_y;i++){
		if(min_dist > dist[i]){
			min_dist = dist[i];
		}
	}
	min_dist = sqrt(min_dist);
	cout<<"reduce minimum distance takes "<<get_time_elapsed(start)<<" ms"<<endl;

	start = get_cur_time();
	cudaFree(d_S);
	cudaFree(d_T);
	cudaFree(d_A);
	cudaFree(d_B);
	cudaFree(d_dist);
	delete dist;
	cout<<"clean spaces takes "<<get_time_elapsed(start)<<" ms"<<endl;


	return min_dist;
}


void SegDist_gpu_fixed_batch(const float *S, const float *T, int batch_size, int batch_num, float *result){

	struct timeval start = get_cur_time();
	float *d_S, *d_T, *d_A, *d_B, *d_dist;

	size_t segment_num = batch_size*batch_num;

	// segment data in device
	cudaMalloc(&d_S, segment_num*6*sizeof(float));
	cudaMalloc(&d_T, segment_num*6*sizeof(float));
	// space for the results in GPU
	cudaMalloc(&d_dist, batch_num*sizeof(float));
	// some temporary space for computation
	cudaMalloc(&d_A, segment_num*3*sizeof(float));
	cudaMalloc(&d_B, segment_num*3*sizeof(float));
	cout<<"allocating space in takes "<<get_time_elapsed(start, true)<<" ms"<<endl;

	cudaMemcpy(d_S, S, segment_num*6*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_T, T, segment_num*6*sizeof(float), cudaMemcpyHostToDevice);
	cout<<"copying data in takes "<<get_time_elapsed(start, true)<<" ms"<<endl;

	// compute the vectors of segments in S and T, save to A and B
	get_vector_kernel<<<batch_num, batch_size>>>(d_S, d_A);
	get_vector_kernel<<<batch_num, batch_size>>>(d_T, d_B);

	cudaDeviceSynchronize();
	cout<<"preprocessing data takes "<<get_time_elapsed(start, true)<<" ms"<<endl;

	// compute the distance in parallel
	SegDist_cuda_new<<<batch_num, batch_size>>>(d_S, d_T, d_A, d_B, d_dist, batch_size);
	cudaDeviceSynchronize();
	cout<<"distances computations takes "<<get_time_elapsed(start, true)<<" ms"<< endl;

	cudaMemcpy(result, d_dist, batch_num*sizeof(float), cudaMemcpyDeviceToHost);
	cout<<"copying result out takes "<<get_time_elapsed(start, true)<<" ms"<<endl;

	cudaFree(d_S);
	cudaFree(d_T);
	cudaFree(d_A);
	cudaFree(d_B);
	cudaFree(d_dist);
	cout<<"clean spaces takes "<<get_time_elapsed(start)<<" ms"<<endl;
}

}
