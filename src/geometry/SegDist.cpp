/*
 * SegDist.cpp
 *
 *  Created on: Nov 23, 2019
 *      Author: teng
 */


#include "geometry.h"


namespace hispeed{

// a helper function to get the distance between seg1
// and seg2. A, B, A_dot_A and B_dot_B are some temporary
// values which are shared by many segment pairs, like A*B and A*C
inline
float SegDist(const float *seg1, const float *seg2,
			  const float *A, const float *B,
			  float A_dot_A, float B_dot_B){

	float X[3], Y[3];

	const float *P = seg1;
	const float *Q = seg2;

	float Tmp[3], A_dot_B, A_dot_T, B_dot_T;

	VmV(Tmp,Q,P);
	A_dot_B = VdotV(A,B);
	A_dot_T = VdotV(A,Tmp);
	B_dot_T = VdotV(B,Tmp);

	if(A_dot_A==0||B_dot_B==0){
		return DBL_MAX;
	}

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
	u = (t*A_dot_B - B_dot_T)/B_dot_B;

	// if u is on segment Q,B, t and u correspond to
	// closest points, otherwise, recompute and
	// clamp t
	if (u <= 0) {
		VcV(Y, Q);
		t = A_dot_T / A_dot_A;
	} else if (u >= 1) {
		VpV(Y, Q, B);
		t = (A_dot_B + A_dot_T) / A_dot_A;
	} else { // on segment
		VpVxS(Y, Q, B, u);
	}

	if (t <= 0) {
		VcV(X, P);
	} else if (t >= 1) {
		VpV(X, P, A);
	} else { // 0<=t<=1
		VpVxS(X, P, A, t);
	}
	VmV(Tmp,X,Y);
	return VdotV(Tmp,Tmp);
}

// check the distance of all size1*size2 pairs
// of segments and return the minimum distance
float SegDist_single(const float *data1, const float *data2,
		size_t size1, size_t size2){
	float local_min = DBL_MAX;
	float *A = new float[size1*3];
	float *B = new float[size2*3];
	float *AdA = new float[size1];
	float *BdB = new float[size2];
	for(int i=0;i<size1;i++){
		VmV(A+i*3, data1+i*6+3, data1+i*6);
		AdA[i] = VdotV(A+i*3, A+i*3);
	}
	for(int i=0;i<size2;i++){
		VmV(B+i*3, data2+i*6+3, data2+i*6);
		BdB[i] = VdotV(B+i*3, B+i*3);
	}

	for(int i=0;i<size1;i++){
		for(int j=0;j<size2;j++){
			const float *cur_S = data1+i*6;
			const float *cur_T = data2+j*6;
			const float *cur_A = A+i*3;
			const float *cur_B = B+j*3;
			float dist = SegDist(cur_S, cur_T, cur_A, cur_B, AdA[i], BdB[j]);
			if(dist < local_min){
				local_min = dist;
			}
		}
	}
	delete A;
	delete B;
	delete AdA;
	delete BdB;
	return local_min;
}

void *SegDist_unit(void *params_void){
	struct dist_param *param = (struct dist_param *)params_void;
	for(int i=0;i<param->batch_num;i++){
		param->dist[i] = SegDist_single(param->data+param->offset_size[4*i]*6,
									    param->data+param->offset_size[4*i+2]*6,
									    param->offset_size[4*i+1],
									    param->offset_size[4*i+3]);
	}
	return NULL;
}

// compute the minimum distance of segment pairs with multiple threads
void SegDist_batch(const float *data, const uint *offset_size, float *result,
				   const uint batch_num, const int num_threads){
	cout<<"starting "<<num_threads<<" threads"<<endl;
	pthread_t threads[num_threads];
	struct dist_param params[num_threads];
	int each_thread = batch_num/num_threads;
	for(int i=0;i<num_threads;i++){
		int start = each_thread*i;
		if(start>=batch_num){
			break;
		}
		params[i].batch_num = min(each_thread, (int)batch_num-start);
		params[i].offset_size = offset_size+start*4;
		params[i].data = data;
		params[i].id = i+1;
		params[i].dist = result+start;
		int rc = pthread_create(&threads[i], NULL, SegDist_unit, (void *)&params[i]);
		if (rc) {
			cout << "Error:unable to create thread," << rc << endl;
			exit(-1);
		}
	}

	for(int i = 0; i < num_threads; i++){
		void *status;
		int rc = pthread_join(threads[i], &status);
		if (rc) {
			cout << "Error:unable to join," << rc << endl;
			exit(-1);
		}
	}
}

}



