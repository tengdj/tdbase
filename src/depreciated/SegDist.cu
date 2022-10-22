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

#include "cuda_util.cuh"
#include "geometry.h"

namespace hispeed{

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

__global__
void SegDist_cuda(const float *data, const uint *offset_size,
				  const float *vec, float *dist, uint cur_offset_1, uint cur_offset_2_start){
	// which batch
	int batch_id = blockIdx.x;
	int cur_offset_2 = threadIdx.x+cur_offset_2_start;

	if(cur_offset_1>=offset_size[batch_id*4+1]){
		return;
	}
	if(cur_offset_2>=offset_size[batch_id*4+3]){
		return;
	}
	uint offset1 = offset_size[batch_id*4];
	uint offset2 = offset_size[batch_id*4+2];

	const float *cur_S = data+6*(offset1+cur_offset_1);
	const float *cur_A = vec+3*(offset1+cur_offset_1);
	const float *cur_T = data+6*(offset2+cur_offset_2);
	const float *cur_B = vec+3*(offset2+cur_offset_2);
	float dd = SegDist_kernel(cur_S, cur_T, cur_A, cur_B);

	if((cur_offset_1==0&&cur_offset_2==0)||dist[batch_id]>dd){
		dist[batch_id]=dd;
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
void get_vector_kernel(float *data, float *vec, int segment_num){
	int id = blockIdx.x*blockDim.x+threadIdx.x;
	if(id<segment_num){
		VmV_d(vec+id*3, data+id*6+3, data+id*6);
	}
}

/*
 * data: contains the segments of the meshes mentioned in this join.
 * offset_size:  contains the offset in the data for each batch, and the sizes of two data sets
 * result: for the returned results for each batch
 * pair_num: number of computed batches
 *
 * */
void SegDist_batch_gpu(gpu_info *gpu, const float *data, const uint *offset_size,
		               float *result, const uint pair_num, const uint segment_num){

	assert(gpu);
	cudaSetDevice(gpu->device_id);
	struct timeval start = get_cur_time();

	// profile the input data
	uint max_size_1 = 0;
	uint max_size_2 = 0;
	for(int i=0;i<pair_num;i++){
		if(offset_size[i*4+1]>max_size_1){
			max_size_1 = offset_size[i*4+1];
		}
		if(offset_size[i*4+3]>max_size_2){
			max_size_2 = offset_size[i*4+3];
		}
	}
	// allocate memory in GPU
	char *cur_d_cuda = gpu->d_data;

	// segment data in device
	float *d_data = (float *)(cur_d_cuda);
	cur_d_cuda += 6*sizeof(float)*segment_num;
	// some temporary space for computation
	float *d_vec = (float *)(cur_d_cuda);
	cur_d_cuda += 3*sizeof(float)*segment_num;
	// space for the results in GPU
	float *d_dist = (float *)(cur_d_cuda);
	cur_d_cuda += sizeof(float)*pair_num;
	// space for the offset and size information in GPU
	uint *d_os = (uint *)(cur_d_cuda);

	CUDA_SAFE_CALL(cudaMemcpy(d_data, data, segment_num*6*sizeof(float), cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_os, offset_size, pair_num*4*sizeof(uint), cudaMemcpyHostToDevice));
	//logt("copying data to GPU", start);

	// compute the vectors of segments in data, save to d_vec
	get_vector_kernel<<<segment_num/1024+1,1024>>>(d_data, d_vec, segment_num);
	check_execution();
	// compute the distance in parallel
	for(uint cur_offset_2=0;cur_offset_2<max_size_2;cur_offset_2+=1024){
		uint dim2 = min(max_size_2-cur_offset_2, 1024);
		for(uint cur_offset_1=0;cur_offset_1<max_size_1;cur_offset_1++){
			SegDist_cuda<<<pair_num, dim2>>>(d_data, d_os, d_vec, d_dist, cur_offset_1, cur_offset_2);
			check_execution();
		}
		//cout<<pair_num<<" "<<cur_offset_2<<" "<<dim2<<endl;
	}

	cudaDeviceSynchronize();
	//logt("distances computations", start);

	CUDA_SAFE_CALL(cudaMemcpy(result, d_dist, pair_num*sizeof(float), cudaMemcpyDeviceToHost));
	//logt("copy data out", start);
}

}
