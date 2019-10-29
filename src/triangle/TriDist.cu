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
inline float SegDist_kernel(const float *seg1, const float *seg2)
{
	float X[3], Y[3];

	const float *P = seg1;
	const float *Q = seg2;

	float A[3], B[3]; // vector of two segments
	VmV(A, seg1+3, seg1);
	VmV(B, seg2+3, seg2);

	float Tmp[3], A_dot_A, B_dot_B, A_dot_B, A_dot_T, B_dot_T;

	VmV(Tmp,Q,P);
	A_dot_A = VdotV(A,A);
	B_dot_B = VdotV(B,B);
	A_dot_B = VdotV(A,B);
	A_dot_T = VdotV(A,Tmp);
	B_dot_T = VdotV(B,Tmp);

	// t parameterizes ray P,A
	// u parameterizes ray Q,B
	float t,u;

	// compute t for the closest point on ray P,A to
	// ray Q,B

	float denom = A_dot_A*B_dot_B - A_dot_B*A_dot_B;
	if(denom == 0){
		t = 0;
	}else{
		t = (A_dot_T*B_dot_B - B_dot_T*A_dot_B) / denom;
	}

	// find u for point on ray Q,B closest to point at t
	if(B_dot_B==0){
		u = 0;
	}else{
		u = (t*A_dot_B - B_dot_T)/B_dot_B;
	}

	// if u is on segment Q,B, t and u correspond to
	// closest points, otherwise, recompute and
	// clamp t
	if (u <= 0) {
		VcV(Y, Q);
		if(A_dot_A==0){
			t = 0;
		}else{
			t = A_dot_T / A_dot_A;
		}

		if (t <= 0) {
			VcV(X, P);
		} else if (t >= 1) {
			VpV(X, P, A);
		} else {
			VpVxS(X, P, A, t);
		}
	} else if (u >= 1) {
		VpV(Y, Q, B);
		if(A_dot_A==0){
			t = 0;
		}else{
			t = (A_dot_B + A_dot_T) / A_dot_A;
		}

		if (t <= 0) {
			VcV(X, P);
		} else if (t >= 1) {
			VpV(X, P, A);
		} else {
			VpVxS(X, P, A, t);
		}
	} else { // on segment
		VpVxS(Y, Q, B, u);
		if (t <= 0) {
			VcV(X, P);
		} else if (t >= 1) {
			VpV(X, P, A);
		} else { // 0<=t<=1
			VpVxS(X, P, A, t);
		}
	}

	VmV(Tmp,X,Y);
	return VdotV(Tmp,Tmp);
}


__global__
void SegDist_cuda(const float *S, const float *T, float *dist){
	int x = blockIdx.x;
	int y = threadIdx.x;
	int id = x*blockDim.x+y;

	float dd = SegDist_kernel(S+x*3, T+y*3);
	if(dist[id]>dd){
		dist[id] = dd;
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
void test_cuda(int *d)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	//int id = blockIdx.x;
	//int dim = blockDim.x;
	d[i] = threadIdx.x;

}

int max_len_x = 200;
int max_len_y = 512;

float SegDist_batch_gpu(const float *S, const float *T, int size1, int size2){
//	float *t_d, *t_d1;
//	float d[10];
//	cudaMalloc(&t_d, 10*sizeof(float));
//	cudaMalloc(&t_d1, 10*sizeof(float));
//	cudaMemset(t_d, 0, 10*sizeof(float));
//
//	get_max<<<2,5>>>(t_d1, t_d);
//	cudaMemcpy(d, t_d, 10*sizeof(int), cudaMemcpyDeviceToHost);
//	for(int i=0;i<10;i++){
//		cout<<i<<" "<<d[i]<<endl;
//	}
//	return 0;

	float *d_x, *d_y, *d_dist;

	int len_x = min(max_len_x, size1);
	int len_y = min(max_len_y, size2);
	float min_dist = DBL_MAX;

	cudaMalloc(&d_x, (size1+1)*3*sizeof(float));
	cudaMalloc(&d_y, (size2+1)*3*sizeof(float));
	cudaMemcpy(d_x, S, (size1+1)*3*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_y, T, (size2+1)*3*sizeof(float), cudaMemcpyHostToDevice);

	cudaMalloc(&d_dist, len_x*len_y*sizeof(float));
	float *dist = new float[len_x*len_y];
	for(int i=0;i<len_x*len_y;i++){
		dist[i] = DBL_MAX;
	}
	cudaMemcpy(d_dist, dist, len_x*len_y*sizeof(float), cudaMemcpyHostToDevice);

	const float *cur_S;
	const float *cur_T;

	cout<<"start computing"<<endl;
	for(int i=0;i<size1;i+=len_x){
		for(int j=0;j<size2;j+=len_y){
			//cout<<i<<" "<<j<<endl;
			cur_S = d_x+i*3;
			cur_T = d_y+j*3;
			int tsize1 = min(len_x, size1-i);
			int tsize2 = min(len_y, size2-j);

			SegDist_cuda<<<tsize1, tsize2>>>(cur_S,cur_T,d_dist);
		}
	}

	cudaMemcpy(dist, d_dist, len_x*len_y*sizeof(float), cudaMemcpyDeviceToHost);
	for(int i=0;i<len_x*len_y;i++){
		if(min_dist > dist[i]){
			min_dist = dist[i];
		}
	}

	cudaFree(d_x);
	cudaFree(d_y);
	cudaFree(d_dist);
	delete dist;

	return sqrt(min_dist);
}
