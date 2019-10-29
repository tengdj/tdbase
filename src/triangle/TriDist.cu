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

//--------------------------------------------------------------------------
// SegPoints()
//
// Returns closest points between an segment pair.
// Implemented from an algorithm described in
//
// Vladimir J. Lumelsky,
// On fast computation of distance between line segments.
// In Information Processing Letters, no. 21, pages 55-61, 1985.
//--------------------------------------------------------------------------

void
SegPoints(float VEC[3],
	  float X[3], float Y[3],             // closest points
	  const float P[3], const float A[3], // seg 1 origin, vector
	  const float Q[3], const float B[3]) // seg 2 origin, vector
{
  float T[3], A_dot_A, B_dot_B, A_dot_B, A_dot_T, B_dot_T;
  float TMP[3];

  VmV(T,Q,P);
  A_dot_A = VdotV(A,A);
  B_dot_B = VdotV(B,B);
  A_dot_B = VdotV(A,B);
  A_dot_T = VdotV(A,T);
  B_dot_T = VdotV(B,T);

  // t parameterizes ray P,A
  // u parameterizes ray Q,B

  float t,u;

  // compute t for the closest point on ray P,A to
  // ray Q,B

  float denom = A_dot_A*B_dot_B - A_dot_B*A_dot_B;

  t = (A_dot_T*B_dot_B - B_dot_T*A_dot_B) / denom;

  // clamp result so t is on the segment P,A

  if ((t < 0) || isnan(t)) t = 0; else if (t > 1) t = 1;

  // find u for point on ray Q,B closest to point at t

  u = (t*A_dot_B - B_dot_T) / B_dot_B;

  // if u is on segment Q,B, t and u correspond to
  // closest points, otherwise, clamp u, recompute and
  // clamp t

  if ((u <= 0) || isnan(u)) {

    VcV(Y, Q);

    t = A_dot_T / A_dot_A;

    if ((t <= 0) || isnan(t)) {
      VcV(X, P);
      VmV(VEC, Q, P);
    }
    else if (t >= 1) {
      VpV(X, P, A);
      VmV(VEC, Q, X);
    }
    else {
      VpVxS(X, P, A, t);
      VcrossV(TMP, T, A);
      VcrossV(VEC, A, TMP);
    }
  }
  else if (u >= 1) {

    VpV(Y, Q, B);

    t = (A_dot_B + A_dot_T) / A_dot_A;

    if ((t <= 0) || isnan(t)) {
      VcV(X, P);
      VmV(VEC, Y, P);
    }
    else if (t >= 1) {
      VpV(X, P, A);
      VmV(VEC, Y, X);
    }
    else {
      VpVxS(X, P, A, t);
      VmV(T, Y, P);
      VcrossV(TMP, T, A);
      VcrossV(VEC, A, TMP);
    }
  }
  else {

    VpVxS(Y, Q, B, u);

    if ((t <= 0) || isnan(t)) {
      VcV(X, P);
      VcrossV(TMP, T, B);
      VcrossV(VEC, B, TMP);
    }
    else if (t >= 1) {
      VpV(X, P, A);
      VmV(T, Q, X);
      VcrossV(TMP, T, B);
      VcrossV(VEC, B, TMP);
    }
    else {
      VpVxS(X, P, A, t);
      VcrossV(VEC, A, B);
      if (VdotV(VEC, T) < 0) {
        VxS(VEC, VEC, -1);
      }
    }
  }
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

float
TriDist(const float *S, const float *T)
{
	return 0;
	bool shown_disjoint = false;

	// closest points
	float P[3];
	float Q[3];

	// some temporary vectors
	float V[3];
	float Z[3];
	// Compute vectors along the 6 sides
	float Sv[3][3], Tv[3][3];

	VmV(Sv[0],S+3,S);
	VmV(Sv[1],S+6,S+3);
	VmV(Sv[2],S,S+6);

	VmV(Tv[0],T+3,T);
	VmV(Tv[1],T+6,T+3);
	VmV(Tv[2],T,T+6);

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
			VmV(V,Q,P);
			float dd = VdotV(V,V);
			if (dd <= mindd){
				mindd = dd;

				// Verify this closest point pair for the segment pairs with minimum distance
				VmV(Z,S+((i+2)%3)*3,P);
				float a = VdotV(Z,VEC);
				VmV(Z,T+((j+2)%3)*3,Q);
				float b = VdotV(Z,VEC);

				// the closest distance of segment pairs is the closest distance of the two triangles
				if ((a <= 0) && (b >= 0)) {
					return sqrt(mindd);
				}

				// otherwise, check the other cases
				// we can use the side product of this calculation
				// to judge whether two triangle joint or not
				float p = VdotV(V, VEC);
				if (a < 0) a = 0;
				if (b > 0) b = 0;
				if ((p - a + b) > 0) shown_disjoint = true;
			}
		}
	}


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

	// First check for case 1

	float Sn[3], Snl;
	VcrossV(Sn,Sv[0],Sv[1]); // Compute normal to S triangle
	Snl = VdotV(Sn,Sn);      // Compute square of length of normal

	// If cross product is long enough,

	if (Snl > 1e-15){
		// Get projection lengths of T points

		float Tp[3];

		VmV(V,S,T);
		Tp[0] = VdotV(V,Sn);

		VmV(V,S,T+3);
		Tp[1] = VdotV(V,Sn);

		VmV(V,S,T+6);
		Tp[2] = VdotV(V,Sn);

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

			VmV(V,T+point*3,S);
			VcrossV(Z,Sn,Sv[0]);
			if (VdotV(V,Z) > 0){
				VmV(V,T+point*3,S+3);
				VcrossV(Z,Sn,Sv[1]);
				if (VdotV(V,Z) > 0) {
					VmV(V,T+point*3,S+6);
					VcrossV(Z,Sn,Sv[2]);
					if (VdotV(V,Z) > 0) {
						// T[point] passed the test - it's a closest point for
						// the T triangle; the other point is on the face of S

						VpVxS(P,T+point*3,Sn,Tp[point]/Snl);
						VcV(Q,T+point*3);
						return sqrt(VdistV2(P,Q));
					}
				}
			}
		}
	}

	float Tn[3], Tnl;
	VcrossV(Tn,Tv[0],Tv[1]);
	Tnl = VdotV(Tn,Tn);

	if (Tnl > 1e-15){

		float Sp[3];
		VmV(V,T,S);
		Sp[0] = VdotV(V,Tn);

		VmV(V,T,S+3);
		Sp[1] = VdotV(V,Tn);

		VmV(V,T,S+6);
		Sp[2] = VdotV(V,Tn);

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

			VmV(V,S+3*point,T);
			VcrossV(Z,Tn,Tv[0]);
			if (VdotV(V,Z) > 0){
				VmV(V,S+3*point,T+3);
				VcrossV(Z,Tn,Tv[1]);
				if (VdotV(V,Z) > 0){
					VmV(V,S+3*point,T+6);
					VcrossV(Z,Tn,Tv[2]);
					if (VdotV(V,Z) > 0){
						VcV(P,S+3*point);
						VpVxS(Q,S+3*point,Tn,Sp[point]/Tnl);
						return sqrt(VdistV2(P,Q));
					}
				}
			}
		}
	}

	// Case 1 can't be shown.
	// If one of these tests showed the triangles disjoint,
	// we assume case 3 or 4, otherwise we conclude case 2,
	// that the triangles overlap.
	if (shown_disjoint){
		return sqrt(mindd);
	}else {
		return 0;
	}
}

/*
 * param dist: head address of space for distances, size = s1*s2*sizeof(float)
 * param P: closest points in S
 * param Q: closest points in T
 * param S: first triangle set size = s1*3*3*sizeof(float)
 * param T: second triangle set size = s2*3*3*sizeof(float)
 * param s1: number of triangles in set S
 * param s2: number of triangles in set T
 *
 * */

struct TriDist_param{
	const float *S;
	const float * T;
	int s1;
	int s2;
	int num_threads;
	int id;
	float dist;
};

void *TriDist_unit(void *params_void){
	struct TriDist_param *param = (struct TriDist_param *)params_void;
	//cout<<"thread "<<param->id<<" is started"<<endl;
	for(int i=0;i<param->s1;i++){
		for(int j=0;j<param->s2;j++){
			const float *cur_S = param->S+i*9;
			const float *cur_T = param->T+j*9;
			float dist = TriDist(cur_S, cur_T);
			if(dist < param->dist){
				param->dist = dist;
			}
		}
	}
	return NULL;
}

float TriDist_batch(const float *S, const float * T, int s1, int s2, int num_threads){
	cout<<"starting "<<num_threads<<" threads"<<endl;

	pthread_t threads[num_threads];
	struct TriDist_param params[num_threads];
	int each_thread = s1/num_threads;
	for(int i=0;i<num_threads;i++){
		int start = each_thread*i;
		if(start>=s1){
			break;
		}
		params[i].s1 = min(start+each_thread, s1)-start;
		params[i].s2 = s2;
		params[i].S = S+start*9;
		params[i].T = T;
		params[i].id = i+1;
		params[i].dist = DBL_MAX;

		int rc = pthread_create(&threads[i], NULL, TriDist_unit, (void *)&params[i]);
		if (rc) {
			cout << "Error:unable to create thread," << rc << endl;
			exit(-1);
		}
	}


	float min = DBL_MAX;
	for(int i = 0; i < num_threads; i++ ){
		void *status;
		int rc = pthread_join(threads[i], &status);
		if (rc) {
			cout << "Error:unable to join," << rc << endl;
			exit(-1);
		}
		if(params[i].dist < min){
			min = params[i].dist;
		}
		//cerr << "Main: completed thread id :" << i <<endl;
	}

	return min;
}

