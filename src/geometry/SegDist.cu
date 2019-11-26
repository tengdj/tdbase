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

#include "geometry.h"

#define CUDA_ERROR_CHECK

#define CudaSafeCall( err ) __cudaSafeCall( err, __FILE__, __LINE__ )
#define CudaCheckError()    __cudaCheckError( __FILE__, __LINE__ )

#define CUDA_EXECUTE(call) \
	do{\
	call; \
	cudaError_t err = cudaGetLastError();\
	if (err != cudaSuccess) printf("Error: %s\n", cudaGetErrorString(err));\
	}while(0);\

#define CUDA_SAFE_CALL(call) \
	do {\
		cudaError_t err = call;\
		if (cudaSuccess != err) {\
			fprintf (stderr, "Cuda error in file '%s' in line %i : %s.\n",\
					__FILE__, __LINE__, cudaGetErrorString(err) );\
			exit(EXIT_FAILURE);\
		}\
	} while (0);

namespace hispeed{

inline void __cudaSafeCall( cudaError err)
{
    if ( cudaSuccess != err )
    {
        fprintf( stderr, "cudaSafeCall() failed: %s\n",
                 cudaGetErrorString( err ) );
        exit( -1 );
    }
    return;
}

inline void check_execution(){
	cudaError_t err = cudaGetLastError();
	if (err != cudaSuccess)
		printf("Error: %s\n", cudaGetErrorString(err));
}

// copy
__device__
inline void
VcV_d(float Vr[3], const float V[3])
{
  Vr[0] = V[0];  Vr[1] = V[1];  Vr[2] = V[2];
}

// minus
__device__
inline void
VmV_d(float Vr[3], const float V1[3], const float V2[3])
{
  Vr[0] = V1[0] - V2[0];
  Vr[1] = V1[1] - V2[1];
  Vr[2] = V1[2] - V2[2];
}

// plus
__device__
inline void
VpV_d(float Vr[3], const float V1[3], const float V2[3])
{
  Vr[0] = V1[0] + V2[0];
  Vr[1] = V1[1] + V2[1];
  Vr[2] = V1[2] + V2[2];
}

// plus after product
__device__
inline void
VpVxS_d(float Vr[3], const float V1[3], const float V2[3], float s)
{
  Vr[0] = V1[0] + V2[0] * s;
  Vr[1] = V1[1] + V2[1] * s;
  Vr[2] = V1[2] + V2[2] * s;
}

// dot product
__device__
inline float
VdotV_d(const float V1[3], const float V2[3])
{
  return (V1[0]*V2[0] + V1[1]*V2[1] + V1[2]*V2[2]);
}

// return the distance of two segments
__device__
inline float SegDist_kernel(const float *S, const float *T,
							const float *A, const float *B)
{

	float t = 0.0, u = 0.0, dist = 0.0, t1 = 0.0;
	float ST[3]; // temporary vector S->T
	VmV_d(ST,T,S);
	float A_dot_A = VdotV_d(A,A);
	float B_dot_B = VdotV_d(B,B);
	if(A_dot_A==0||B_dot_B==0){
		return DBL_MAX;
	}
	float A_dot_B = VdotV_d(A,B);
	float A_dot_ST = VdotV_d(A,ST);
	float B_dot_ST = VdotV_d(B,ST);

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
void SegDist_cuda_new(const float *data, const long *offset_size,
				  	  const float *vec, float *dist){
	// computing which segment in set1 of voxel with
	// voxel_id
	int voxel_id = blockIdx.x;
	int segment_id = threadIdx.x;
	long offset = offset_size[voxel_id*3];
	long size1 = offset_size[voxel_id*3+1];
	long size2 = offset_size[voxel_id*3+2];

	// update the pointers for current thread
	const float *cur_S = data+6*(offset+segment_id);
	const float *cur_T = data+6*(offset+size1);
	const float *cur_A = vec+3*(offset+segment_id);
	const float *cur_B = vec+3*(offset+size1);

	float min_dist = DBL_MAX;
	// go over all the segments in set2
	for(int i = 0;i<size2;i++){
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
	for(int i=1;i<size1;i++){
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
void get_vector_kernel(float *data, float *vec){
	int id = threadIdx.x;
	VmV_d(vec+id*3, data+id*6+3, data+id*6);
}

/*
 * data: contains the data of the objects, pact set1 and set2 together.
 * offset_size:  contains the offset in the data for each batch, and the sizes of two data sets
 * result: for the returned results for each batch
 * batch_num: number of computed batches
 *
 * */
void SegDist_batch_gpu(const float *data, const long *offset_size, float *result, const int batch_num){

	size_t segment_num = 0;
	// we use the maximum size of elements
	// in set 1 as one computing dimension
	long batch_size = LONG_MIN;
	for(int i=0;i<batch_num;i++){
		segment_num += offset_size[3*i+1];
		segment_num += offset_size[3*i+2];
		if(batch_size<offset_size[3*i+1]){
			batch_size = offset_size[3*i+1];
		}
	}
	struct timeval start = get_cur_time();
	float *d_data, *d_vec, *d_dist;
	long *d_os;

	// segment data in device
	CUDA_SAFE_CALL(cudaMalloc(&d_data, segment_num*6*sizeof(float)));
	// space for the offset and size information in GPU
	CUDA_SAFE_CALL(cudaMalloc(&d_os, 3*batch_num*sizeof(long)));
	// space for the results in GPU
	CUDA_SAFE_CALL(cudaMalloc(&d_dist, batch_num*sizeof(float)));
	// some temporary space for computation
	CUDA_SAFE_CALL(cudaMalloc(&d_vec, segment_num*3*sizeof(float)));
	report_time("allocating space in GPU", start);

	CUDA_SAFE_CALL(cudaMemcpy(d_data, data, segment_num*6*sizeof(float), cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_os, offset_size, batch_num*3*sizeof(long), cudaMemcpyHostToDevice));
	report_time("copying data to GPU", start);
	cudaError_t err = cudaGetLastError();
	if (err != cudaSuccess) printf("Error: %s\n", cudaGetErrorString(err));
	// compute the vectors of segments in data, save to d_vec
	get_vector_kernel<<<1,segment_num>>>(d_data, d_vec);
	check_execution();
	int num_per_iter = 400;
	int cur_iter = 0;
	do{
		cur_iter++;
		if(cur_iter%num_per_iter==0||cur_iter==batch_num){
			int step = (cur_iter==batch_num)?(cur_iter%num_per_iter):num_per_iter;
			long cur_offset = offset_size[(cur_iter-step)*3];
			float *cur_data = d_data+cur_offset*6;
			float *cur_vec = d_vec+cur_offset*3;
			long *cur_os = d_os+(cur_iter-step)*3;
			float *cur_dist = d_dist+(cur_iter-step);
			SegDist_cuda_new<<<step, batch_size>>>(cur_data, cur_os, cur_vec, cur_dist);
			check_execution();
		}
		break;
	}while(cur_iter<batch_num);
	// compute the distance in parallel
	cudaDeviceSynchronize();
	report_time("distances computations", start);

	CUDA_SAFE_CALL(cudaMemcpy(result, d_dist, batch_num*sizeof(float), cudaMemcpyDeviceToHost));
	CUDA_SAFE_CALL(cudaFree(d_data));
	CUDA_SAFE_CALL(cudaFree(d_vec));
	CUDA_SAFE_CALL(cudaFree(d_os));
	CUDA_SAFE_CALL(cudaFree(d_dist));
	report_time("copy data out and clean", start);

}

}
