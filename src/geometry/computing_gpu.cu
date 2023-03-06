/*
 * TriDist.cu
 *
 *  Created on: Oct 22, 2022
 *      Author: teng
 */

#include <pthread.h>
#include "geometry.h"
#include "cuda_util.cuh"
namespace hispeed{


/*
 *
 * get the closest points between segments
 *
 * */
__device__
inline void
SegPoints(float VEC[3],
	  float X[3], float Y[3],             // closest points
	  const float P[3], const float A[3], // seg 1 origin, vector
	  const float Q[3], const float B[3]) // seg 2 origin, vector
{
  float T[3], A_dot_A, B_dot_B, A_dot_B, A_dot_T, B_dot_T;
  float TMP[3];

  VmV_d(T,Q,P);
  A_dot_A = VdotV_d(A,A);
  B_dot_B = VdotV_d(B,B);
  A_dot_B = VdotV_d(A,B);
  A_dot_T = VdotV_d(A,T);
  B_dot_T = VdotV_d(B,T);
  assert(A_dot_A!=0&&B_dot_B!=0);

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
    VcV_d(Y, Q);
    if(A_dot_A==0){
    	t = 0;
    }else{
        t = A_dot_T / A_dot_A;
    }

    if (t <= 0) {
      VcV_d(X, P);
      VmV_d(VEC, Q, P);
    }
    else if (t >= 1) {
      VpV_d(X, P, A);
      VmV_d(VEC, Q, X);
    }
    else {
      VpVxS_d(X, P, A, t);
      VcrossV_d(TMP, T, A);
      VcrossV_d(VEC, A, TMP);
    }
  }
  else if (u >= 1) {
    VpV_d(Y, Q, B);
    if(A_dot_A==0){
    	t = 0;
    }else{
        t = (A_dot_B + A_dot_T) / A_dot_A;
    }

    if (t <= 0) {
      VcV_d(X, P);
      VmV_d(VEC, Y, P);
    }
    else if (t >= 1) {
      VpV_d(X, P, A);
      VmV_d(VEC, Y, X);
    }
    else {
      VpVxS_d(X, P, A, t);
      VmV_d(T, Y, P);
      VcrossV_d(TMP, T, A);
      VcrossV_d(VEC, A, TMP);
    }
  }
  else { // on segment

    VpVxS_d(Y, Q, B, u);

    if (t <= 0) {
      VcV_d(X, P);
      VcrossV_d(TMP, T, B);
      VcrossV_d(VEC, B, TMP);
    }
    else if (t >= 1) {
      VpV_d(X, P, A);
      VmV_d(T, Q, X);
      VcrossV_d(TMP, T, B);
      VcrossV_d(VEC, B, TMP);
    }
    else { // 0<=t<=1
      VpVxS_d(X, P, A, t);
      VcrossV_d(VEC, A, B);
      if (VdotV_d(VEC, T) < 0) {
        VxS_d(VEC, VEC, -1);
      }
    }
  }
}

/*
 * check the segments of a triangle to see if the closest
 * points can be found on a segment pair, which covers
 * almost all cases
 * */
__device__
inline float
TriDist_seg(const float *S, const float *T,
		bool &shown_disjoint, bool &closest_find){

	// closest points
	float P[3];
	float Q[3];

	// some temporary vectors
	float V[3];
	float Z[3];
	// Compute vectors along the 6 sides
	float Sv[3][3], Tv[3][3];

	VmV_d(Sv[0],S+3,S);
	VmV_d(Sv[1],S+6,S+3);
	VmV_d(Sv[2],S,S+6);

	VmV_d(Tv[0],T+3,T);
	VmV_d(Tv[1],T+6,T+3);
	VmV_d(Tv[2],T,T+6);

	// For each edge pair, the vector connecting the closest points
	// of the edges defines a slab (parallel planes at head and tail
	// enclose the slab). If we can show that the off-edge vertex of
	// each triangle is outside of the slab, then the closest points
	// of the edges are the closest points for the triangles.
	// Even if these tests fail, it may be helpful to know the closest
	// points found, and whether the triangles were shown disjoint

	float mindd = DBL_MAX; // Set first minimum safely high
	float VEC[3];
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			// Find closest points on edges i & j, plus the
			// vector (and distance squared) between these points
			SegPoints(VEC,P,Q,S+i*3,Sv[i],T+j*3,Tv[j]);
			VmV_d(V,Q,P);
			float dd = VdotV_d(V,V);
			if (dd <= mindd){
				mindd = dd;

				// Verify this closest point pair for the segment pairs with minimum distance
				VmV_d(Z,S+((i+2)%3)*3,P);
				float a = VdotV_d(Z,VEC);
				VmV_d(Z,T+((j+2)%3)*3,Q);
				float b = VdotV_d(Z,VEC);

				// the closest distance of segment pairs is the closest distance of the two triangles
				if ((a <= 0) && (b >= 0)) {
					closest_find = true;
					return sqrt(mindd);
				}

				// otherwise, check the other cases
				// we can use the side product of this calculation
				// to judge whether two triangle joint or not
				float p = VdotV_d(V, VEC);
				if (a < 0) a = 0;
				if (b > 0) b = 0;
				if ((p - a + b) > 0) shown_disjoint = true;
			}
		}
	}

	// not sure is the case
	closest_find = false;
	return mindd;
}

/*
 * any other cases for triangle distance calculation besides the
 * closest points reside on segments
 *
 * */
__device__
inline float
TriDist_other(const float *S, const float *T, bool &shown_disjoint)
{

	// closest points
	float P[3];
	float Q[3];

	// some temporary vectors
	float V[3];
	float Z[3];
	// Compute vectors along the 6 sides
	float Sv[3][3], Tv[3][3];

	VmV_d(Sv[0],S+3,S);
	VmV_d(Sv[1],S+6,S+3);
	VmV_d(Sv[2],S,S+6);

	VmV_d(Tv[0],T+3,T);
	VmV_d(Tv[1],T+6,T+3);
	VmV_d(Tv[2],T,T+6);

	// First check for case 1

	float Sn[3], Snl;
	VcrossV_d(Sn,Sv[0],Sv[1]); // Compute normal to S triangle
	Snl = VdotV_d(Sn,Sn);      // Compute square of length of normal

	// If cross product is long enough,

	if (Snl > 1e-15){
		// Get projection lengths of T points

		float Tp[3];

		VmV_d(V,S,T);
		Tp[0] = VdotV_d(V,Sn);

		VmV_d(V,S,T+3);
		Tp[1] = VdotV_d(V,Sn);

		VmV_d(V,S,T+6);
		Tp[2] = VdotV_d(V,Sn);

		// If Sn is a separating direction,
		// find point with smallest projection

		int point = -1;
		if ((Tp[0] > 0) && (Tp[1] > 0) && (Tp[2] > 0)) {
			if (Tp[0] < Tp[1]) point = 0; else point = 1;
			if (Tp[2] < Tp[point]) point = 2;
		} else if ((Tp[0] < 0) && (Tp[1] < 0) && (Tp[2] < 0)) {
			if (Tp[0] > Tp[1]) point = 0; else point = 1;
			if (Tp[2] > Tp[point]) point = 2;
		}

		// If Sn is a separating direction,

		if (point >= 0){
			shown_disjoint = true;
			// Test whether the point found, when projected onto the
			// other triangle, lies within the face.

			VmV_d(V,T+point*3,S);
			VcrossV_d(Z,Sn,Sv[0]);
			if (VdotV_d(V,Z) > 0){
				VmV_d(V,T+point*3,S+3);
				VcrossV_d(Z,Sn,Sv[1]);
				if (VdotV_d(V,Z) > 0) {
					VmV_d(V,T+point*3,S+6);
					VcrossV_d(Z,Sn,Sv[2]);
					if (VdotV_d(V,Z) > 0) {
						// T[point] passed the test - it's a closest point for
						// the T triangle; the other point is on the face of S

						VpVxS_d(P,T+point*3,Sn,Tp[point]/Snl);
						VcV_d(Q,T+point*3);
						return sqrt(VdistV2_d(P,Q));
					}
				}
			}
		}
	}

	float Tn[3], Tnl;
	VcrossV_d(Tn,Tv[0],Tv[1]);
	Tnl = VdotV_d(Tn,Tn);

	if (Tnl > 1e-15){

		float Sp[3];
		VmV_d(V,T,S);
		Sp[0] = VdotV_d(V,Tn);

		VmV_d(V,T,S+3);
		Sp[1] = VdotV_d(V,Tn);

		VmV_d(V,T,S+6);
		Sp[2] = VdotV_d(V,Tn);

		int point = -1;
		if ((Sp[0] > 0) &&
				(Sp[1] > 0) && (Sp[2] > 0)) {
			if (Sp[0] < Sp[1]) point = 0; else point = 1;
			if (Sp[2] < Sp[point]) point = 2;
		} else if ((Sp[0] < 0) &&
				(Sp[1] < 0) && (Sp[2] < 0)) {
			if (Sp[0] > Sp[1]) point = 0; else point = 1;
			if (Sp[2] > Sp[point]) point = 2;
		}

		if (point >= 0){
			shown_disjoint = true;

			VmV_d(V,S+3*point,T);
			VcrossV_d(Z,Tn,Tv[0]);
			if (VdotV_d(V,Z) > 0){
				VmV_d(V,S+3*point,T+3);
				VcrossV_d(Z,Tn,Tv[1]);
				if (VdotV_d(V,Z) > 0){
					VmV_d(V,S+3*point,T+6);
					VcrossV_d(Z,Tn,Tv[2]);
					if (VdotV_d(V,Z) > 0){
						VcV_d(P,S+3*point);
						VpVxS_d(Q,S+3*point,Tn,Sp[point]/Tnl);
						return sqrt(VdistV2_d(P,Q));
					}
				}
			}
		}
	}

	// not the case
	return -1;
}

//--------------------------------------------------------------------------
// TriDist()
//
// Computes the closest points on two triangles, and returns the
// distance between them.
//
// S and T are the triangles, stored tri[point][dimension].
//
// If the triangles are disjoint, P and Q give the closest points of
// S and T respectively. However, if the triangles overlap, P and Q
// are basically a random pair of points from the triangles, not
// coincident points on the intersection of the triangles, as might
// be expected.
//--------------------------------------------------------------------------

__device__
float
TriDist_kernel(const float *S, const float *T)
{
	bool shown_disjoint = false;
	bool closest_find = false;
	float mindd_seg = TriDist_seg(S, T, shown_disjoint, closest_find);
	if(closest_find){// the closest points are one segments, simply return
		return mindd_seg;
	}else{
		// No edge pairs contained the closest points.
		// either:
		// 1. one of the closest points is a vertex, and the
		//    other point is interior to a face.
		// 2. the triangles are overlapping.
		// 3. an edge of one triangle is parallel to the other's face. If
		//    cases 1 and 2 are not true, then the closest points from the 9
		//    edge pairs checks above can be taken as closest points for the
		//    triangles.
		// 4. possibly, the triangles were degenerate.  When the
		//    triangle points are nearly colinear or coincident, one
		//    of above tests might fail even though the edges tested
		//    contain the closest points.
		float mindd_other = TriDist_other(S, T, shown_disjoint);
		if(mindd_other != -1){ // is the case
			return mindd_other;
		}
	}

	// Case 1 can't be shown.
	// If one of these tests showed the triangles disjoint,
	// we assume case 3 or 4, otherwise we conclude case 2,
	// that the triangles overlap.
	if (shown_disjoint){
		return sqrt(mindd_seg);
	}else {
		return 0;
	}
}

__global__
void TriDist_cuda(const float *data, const uint *offset_size, result_container *dist, uint cur_offset_1, uint cur_offset_2_start){
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

	const float *cur_S = data+9*(offset1+cur_offset_1);
	const float *cur_T = data+9*(offset2+cur_offset_2);
	float dd = TriDist_kernel(cur_S, cur_T);

	if((cur_offset_1==0&&cur_offset_2==0)||dist[batch_id].result.distance>dd){
		dist[batch_id].result.distance = dd;
		dist[batch_id].p1 = cur_offset_1;
		dist[batch_id].p2 = cur_offset_2;
	}
}

__global__
void TriInt_cuda(const float *data, const uint *offset_size, result_container *intersect, uint cur_offset_1, uint cur_offset_2_start){

	int batch_id = blockIdx.x;
	uint cur_offset_2 = threadIdx.x+cur_offset_2_start;

	if(cur_offset_1>=offset_size[batch_id*4+1]){
		return;
	}
	if(cur_offset_2>=offset_size[batch_id*4+3]){
		return;
	}

	// determined
	if(intersect[batch_id].result.intersected){
		return;
	}

	uint offset1 = offset_size[batch_id*4];
	uint offset2 = offset_size[batch_id*4+2];
	const float *cur_S = data+9*(offset1+cur_offset_1);
	const float *cur_T = data+9*(offset2+cur_offset_2);

	float dd = TriDist_kernel(cur_S, cur_T);
	if(dd==0.0){
		intersect[batch_id].result.intersected = 1;
		intersect[batch_id].p1 = cur_offset_1;
		intersect[batch_id].p2 = cur_offset_2;
		return;
	}
}


/*
 * data: contains the triangles of the meshes in this join.
 * offset_size:  contains the offset in the data for each batch, and the sizes of two data sets
 * result: for the returned results for each batch
 * pair_num: number of computed batches
 *
 * */
void TriDist_batch_gpu(gpu_info *gpu, const float *data, const uint *offset_size,
		               result_container *result, const uint pair_num, const uint element_num){

	assert(gpu);
	cudaSetDevice(gpu->device_id);
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
	cur_d_cuda += 9*sizeof(float)*element_num;
	// space for the results in GPU
	result_container *d_dist = (result_container *)(cur_d_cuda);
	cur_d_cuda += sizeof(result_container)*pair_num;
	// space for the offset and size information in GPU
	uint *d_os = (uint *)(cur_d_cuda);

	CUDA_SAFE_CALL(cudaMemcpy(d_data, data, element_num*9*sizeof(float), cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_os, offset_size, pair_num*4*sizeof(uint), cudaMemcpyHostToDevice));
	//logt("copying data to GPU", start);

	// compute the distance in parallel
	for(uint cur_offset_2=0;cur_offset_2<max_size_2;cur_offset_2+=1024){
		uint dim2 = min(max_size_2-cur_offset_2, (uint)512);
		for(uint cur_offset_1=0;cur_offset_1<max_size_1;cur_offset_1++){
			TriDist_cuda<<<pair_num, dim2>>>(d_data, d_os, d_dist, cur_offset_1, cur_offset_2);
			check_execution();
		}
		//cout<<pair_num<<" "<<cur_offset_2<<" "<<dim2<<endl;
	}
	cudaDeviceSynchronize();
	//logt("distances computations", start);

	CUDA_SAFE_CALL(cudaMemcpy(result, d_dist, pair_num*sizeof(result_container), cudaMemcpyDeviceToHost));
	//logt("copy data out", start);
}

__global__
void clear_resultset(result_container *result, uint pairnum){

	int id = blockIdx.x*blockDim.x+threadIdx.x;
	if(id>=pairnum){
		return;
	}
	result[id].result.intersected = 0;
}

void TriInt_batch_gpu(gpu_info *gpu, const float *data, const uint *offset_size,
		result_container *intersection, const uint pair_num, const uint triangle_num){

	assert(gpu);
	cudaSetDevice(gpu->device_id);
	struct timeval start = get_cur_time();
	// allocate memory in GPU
	char *cur_d_cuda = gpu->d_data;

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

	//log("%ld %ld", max_size_1, max_size_2);
	// segment data in device
	float *d_data = (float *)(cur_d_cuda);
	cur_d_cuda += 9*triangle_num*sizeof(float);
	// space for the results in GPU
	result_container *d_intersect = (result_container *)(cur_d_cuda);
	cur_d_cuda += sizeof(result_container)*pair_num;
	// space for the offset and size information in GPU
	uint *d_os = (uint *)(cur_d_cuda);

	CUDA_SAFE_CALL(cudaMemcpy(d_data, data, triangle_num*9*sizeof(float), cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_os, offset_size, pair_num*4*sizeof(uint), cudaMemcpyHostToDevice));
	//logt("copying data to GPU", start);

	clear_resultset<<<pair_num/1024+1, 1024>>>(d_intersect, pair_num);

	// check the intersection
	for(uint cur_offset_2=0;cur_offset_2<max_size_2;cur_offset_2+=1024){
		uint dim2 = min(max_size_2-cur_offset_2, (uint)1024);
		for(uint cur_offset_1=0;cur_offset_1<max_size_1;cur_offset_1++){
			TriInt_cuda<<<pair_num, dim2>>>(d_data, d_os, d_intersect, cur_offset_1,cur_offset_2);
			check_execution();
		}
	}
	check_execution();

	//cout<<pair_num<<" "<<triangle_num<<" "<<sizeof(uint)<<endl;
	cudaDeviceSynchronize();
	//logt("distances computations", start);

	CUDA_SAFE_CALL(cudaMemcpy(intersection, d_intersect, pair_num*sizeof(result_container), cudaMemcpyDeviceToHost));
	//logt("copy data out", start);
}


}

