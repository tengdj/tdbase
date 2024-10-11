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

#include <pthread.h>
#include "geometry.h"

namespace tdbase{


/*
 *
 * get the closest points between segments
 *
 * */
inline void
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
    VcV(Y, Q);
    if(A_dot_A==0){
    	t = 0;
    }else{
        t = A_dot_T / A_dot_A;
    }

    if (t <= 0) {
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
    if(A_dot_A==0){
    	t = 0;
    }else{
        t = (A_dot_B + A_dot_T) / A_dot_A;
    }

    if (t <= 0) {
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
  else { // on segment

    VpVxS(Y, Q, B, u);

    if (t <= 0) {
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
    else { // 0<=t<=1
      VpVxS(X, P, A, t);
      VcrossV(VEC, A, B);
      if (VdotV(VEC, T) < 0) {
        VxS(VEC, VEC, -1);
      }
    }
  }
}

/*
 * check the segments of a triangle to see if the closest
 * points can be found on a segment pair, which covers
 * almost all cases
 * */
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
					closest_find = true;
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

	// not sure is the case
	closest_find = false;
	return mindd;
}

/*
 * any other cases for triangle distance calculation besides the
 * closest points reside on segments
 *
 * */
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

	VmV(Sv[0],S+3,S);
	VmV(Sv[1],S+6,S+3);
	VmV(Sv[2],S,S+6);

	VmV(Tv[0],T+3,T);
	VmV(Tv[1],T+6,T+3);
	VmV(Tv[2],T,T+6);

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

float TriDist(const float *S, const float *T)
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

result_container MeshDist(const float *data1, const float *data2, size_t size1, size_t size2, const float *hausdorff1, const float *hausdorff2){
	result_container ret;
	ret.distance = DBL_MAX;
	ret.min_dist = DBL_MAX;
	ret.max_dist = DBL_MAX;
	for(size_t i=0;i<size1;i++){
		for(size_t j=0;j<size2;j++){
			// get distance of current triangle pair
			// each triangle contains three points
			float dist = TriDist(data1+i*9, data2+j*9);
			if(dist < ret.distance){
				ret.distance = dist;
				ret.p1 = i;
				ret.p2 = j;
			}

			if(hausdorff1 == NULL || hausdorff2 == NULL){
				continue;
			}

			// with hausdorff distances under consideration
			float phdist1 = *(hausdorff1+2*i);
			float phdist2 = *(hausdorff2+2*j);
			float hdist1 = *(hausdorff1+2*i+1);
			float hdist2 = *(hausdorff2+2*j+1);

			float low_dist = std::max(dist-phdist1-phdist2, (float)0.0);
			float high_dist = dist+hdist1+hdist2;
			ret.min_dist = min(ret.min_dist, low_dist);
			ret.max_dist = min(ret.max_dist, high_dist);

//			const float *tdata1 = data1+i*9;
//			const float *tdata2 = data2+j*9;
//			printf("(%f %f %f, %f %f %f, %f %f %f) (%f %f %f, %f %f %f, %f %f %f) %f\n"
//					,*tdata1,*(tdata1+1),*(tdata1+2),*(tdata1+3),*(tdata1+4),*(tdata1+5),*(tdata1+6),*(tdata1+7),*(tdata1+8)
//					,*tdata2,*(tdata2+1),*(tdata2+2),*(tdata2+3),*(tdata2+4),*(tdata2+5),*(tdata2+6),*(tdata2+7),*(tdata2+8)
//					,dist);
		}
	}
	return ret;
}

result_container MeshInt(const float *data1, const float *data2, size_t size1, size_t size2, const float *hausdorff1, const float *hausdorff2){
	result_container res;
	res.intersected = false;

	res.distance = DBL_MAX;

	if(!hausdorff1 || !hausdorff2){
		for(size_t i=0;i<size1;i++){
			for(size_t j=0;j<size2;j++){
				if(TriInt(data1+9*i, data2+9*j)){
					res.intersected = true;
					res.p1 = i;
					res.p2 = j;
					return res;
				}
			}
		}
		return res;
	}

	for(size_t i=0;i<size1;i++){
		for(size_t j=0;j<size2;j++){
			float dist = TriDist(data1+9*i, data2+9*j);
			float phdist1 = 0;
			float phdist2 = 0;
			if(hausdorff1 != NULL && hausdorff2 != NULL){
				phdist1 = *(hausdorff1+2*i);
				phdist2 = *(hausdorff2+2*j);
			}
			res.distance = min(res.distance, dist - phdist1 - phdist2);

			if(dist==0) {
				res.intersected = true;
				res.p1 = i;
				res.p2 = j;
				return res;
			}
		}
	}

	return res;
}

float PointTriangleDist(const float *point, const float *triangle)
{
	// The member result.sqrDistance is set in each block of the
	// nested if-then-else statements. The remaining members are all
	// set at the end of the function.

	float const zero = static_cast<float>(0);
	float const one = static_cast<float>(1);
	float const two = static_cast<float>(2);

	float diff[3];
	float edge0[3];
	float edge1[3];
	const float *trianglev0 = triangle;
	const float *trianglev1 = triangle+3;
	const float *trianglev2 = triangle+6;
	VmV(diff, trianglev0, point);
	VmV(edge0, trianglev1, trianglev0);
	VmV(edge1, trianglev2, trianglev0);

	float a00 = VdotV(edge0, edge0);
	float a01 = VdotV(edge0, edge1);
	float a11 = VdotV(edge1, edge1);
	float b0 = VdotV(diff, edge0);
	float b1 = VdotV(diff, edge1);
	float det = std::max(a00 * a11 - a01 * a01, (float)0.0);
	float s = a01 * b1 - a11 * b0;
	float t = a01 * b0 - a00 * b1;

	// some bad triangles
	if(a00==0.0||a01==0.0||a11==0.0||det==0.0){
		return sqrt(VdotV(diff, diff));
	}

	if(a00==0.0||a01==0.0||a11==0.0||det==0.0){
		printf("%f %f %f, %f %f %f, %f %f %f\n",
				*(triangle+0)
				,*(triangle+1)
				,*(triangle+2)
				,*(triangle+3)
				,*(triangle+4)
				,*(triangle+5)
				,*(triangle+6)
				,*(triangle+7)
				,*(triangle+8));
		printf("%f %f %f %f\n", a00, a01, a11, a00 * a11 - a01 * a01);
	}

	if (s + t <= det)
	{
		if (s < 0.0)
		{
			if (t < 0.0)  // region 4
			{
				if (b0 < 0.0)
				{
					t = zero;
					if (-b0 >= a00)
					{
						s = one;
					}
					else
					{
						s = -b0 / a00;
					}
				}
				else
				{
					s = zero;
					if (b1 >= zero)
					{
						t = zero;
					}
					else if (-b1 >= a11)
					{
						t = one;
					}
					else
					{
						t = -b1 / a11;
					}
				}
			}
			else  // region 3
			{
				s = zero;
				if (b1 >= zero)
				{
					t = zero;
				}
				else if (-b1 >= a11)
				{
					t = one;
				}
				else
				{
					t = -b1 / a11;
				}
			}
		}
		else if (t < zero)  // region 5
		{
			t = zero;
			if (b0 >= zero)
			{
				s = zero;
			}
			else if (-b0 >= a00)
			{
				s = one;
			}
			else
			{
				s = -b0 / a00;
			}
		}
		else  // region 0
		{
			// minimum at interior point
			s /= det;
			t /= det;
		}
	}
	else
	{
		float tmp0{}, tmp1{}, numer{}, denom{};

		if (s < zero)  // region 2
		{
			tmp0 = a01 + b0;
			tmp1 = a11 + b1;
			if (tmp1 > tmp0)
			{
				numer = tmp1 - tmp0;
				denom = a00 - two * a01 + a11;
				assert(denom!=0.0);
				if (numer >= denom)
				{
					s = one;
					t = zero;
				}
				else
				{
					assert(denom>0);
					s = numer / denom;
					t = one - s;
				}
			}
			else
			{
				s = zero;
				if (tmp1 <= zero)
				{
					t = one;
				}
				else if (b1 >= zero)
				{
					t = zero;
				}
				else
				{
					t = -b1 / a11;
				}
			}
		}
		else if (t < zero)  // region 6
		{
			tmp0 = a01 + b1;
			tmp1 = a00 + b0;
			if (tmp1 > tmp0)
			{
				numer = tmp1 - tmp0;
				denom = a00 - two * a01 + a11;
				if (numer >= denom)
				{
					t = one;
					s = zero;
				}
				else
				{
					assert(denom>0);
					t = numer / denom;
					s = one - t;
				}
			}
			else
			{
				t = zero;
				if (tmp1 <= zero)
				{
					s = one;
				}
				else if (b0 >= zero)
				{
					s = zero;
				}
				else
				{
					s = -b0 / a00;
				}
			}
		}
		else  // region 1
		{
			numer = a11 + b1 - a01 - b0;
			if (numer <= zero)
			{
				s = zero;
				t = one;
			}
			else
			{
				denom = a00 - two * a01 + a11;
				if (numer >= denom)
				{
					s = one;
					t = zero;
				}
				else
				{
					assert(denom!=0.0);
					s = numer / denom;
					t = one - s;
				}
			}
		}
	}

	float closest[3];
	VpVxS(closest, trianglev0, edge0, s);
	VpVxS(closest, closest, edge1, t);
	VmV(closest, closest, point);
	return sqrt(VdotV(closest, closest));
}

void compute_normal(float *Norm, const float *triangle){
	const float *trianglev0 = triangle;
	const float *trianglev1 = triangle+3;
	const float *trianglev2 = triangle+6;
	float A[3];
	float B[3];
	VmV(A, trianglev1, trianglev0);
	VmV(B, trianglev2, trianglev0);
	VcrossV(Norm, A, B);
	VdS(Norm, Norm, sqrt(VdotV(Norm, Norm)));
}

void project_points_to_triangle_plane(const float *point, const float *triangle, float projected_point[3]){
	const float *trianglev0 = triangle;
	const float *trianglev1 = triangle+3;
	const float *trianglev2 = triangle+6;
	float A[3];
	float B[3];
	VmV(A, trianglev1, trianglev0);
	VmV(B, trianglev2, trianglev0);
	float Norm[3];
	VcrossV(Norm, A, B);
	VdS(Norm, Norm, sqrt(VdotV(Norm, Norm)));
	float d = -VdotV(Norm, trianglev0);

//	cout<<Norm[0]*trianglev0[0]+Norm[1]*trianglev0[1]+Norm[2]*trianglev0[2]+d<<endl;
//	cout<<Norm[0]*trianglev1[0]+Norm[1]*trianglev1[1]+Norm[2]*trianglev1[2]+d<<endl;
//	cout<<Norm[0]*trianglev2[0]+Norm[1]*trianglev2[1]+Norm[2]*trianglev2[2]+d<<endl;
	VmV(A, point, trianglev0);
	float dist = VdotV(A, Norm);
	VpVxS(projected_point, point, Norm, -dist);
}

bool PointInTriangleCylinder(const float *point, const float *triangle)
{
	const float *trianglev0 = triangle;
	const float *trianglev1 = triangle+3;
	const float *trianglev2 = triangle+6;
	float A[3];
	float B[3];
	VmV(A, trianglev1, trianglev0);
	VmV(B, trianglev2, trianglev0);
	float Norm[3];
	VcrossV(Norm, A, B);
	VdS(Norm, Norm, sqrt(VdotV(Norm, Norm)));
	float d = -VdotV(Norm, trianglev0);

//	cout<<Norm[0]*trianglev0[0]+Norm[1]*trianglev0[1]+Norm[2]*trianglev0[2]+d<<endl;
//	cout<<Norm[0]*trianglev1[0]+Norm[1]*trianglev1[1]+Norm[2]*trianglev1[2]+d<<endl;
//	cout<<Norm[0]*trianglev2[0]+Norm[1]*trianglev2[1]+Norm[2]*trianglev2[2]+d<<endl;
	VmV(A, point, trianglev0);
	float dist = VdotV(A, Norm);

	float projected_point[3];
	VpVxS(projected_point, point, Norm, -dist);
	//cout<<projected_point[0]<<" "<<projected_point[1]<<" "<<projected_point[2]<<endl;

	float v0v1[3];
	float v0v2[3];
	float v1v2[3];
	float v0vp[3];
	float v1vp[3];
	VmV(v0v1,trianglev1,trianglev0);
	VmV(v0v2,trianglev2,trianglev0);
	VmV(v1v2,trianglev2,trianglev1);
	VmV(v0vp,projected_point,trianglev0);
	VmV(v1vp,projected_point,trianglev1);

	VcrossV(A, v0v1, v0v2);
	float a = sqrt(VdotV(A,A));

	VcrossV(A, v0v1, v0vp);
	float a1 = sqrt(VdotV(A,A));

	VcrossV(A, v0v2, v0vp);
	float a2 = sqrt(VdotV(A,A));

	VcrossV(A, v1v2, v1vp);
	float a3 = sqrt(VdotV(A,A));

	//cout<<a1+a2+a3<<" "<<a<<endl;
	return abs(a1+a2+a3 - a)<a/10000000.0;

	VmV(A,trianglev1,trianglev0);
	VmV(B,projected_point,trianglev0);
	cout<<"1: "<<VdotV(A,B)<<endl;
	if(VdotV(A,B)<=0.0){
		return false;
	}
	VmV(A,trianglev2,trianglev1);
	VmV(B,projected_point,trianglev1);
	cout<<"2: "<<VdotV(A,B)<<endl;
	if(VdotV(A,B)<=0.0){
		return false;
	}
	VmV(A,trianglev0,trianglev2);
	VmV(B,projected_point,trianglev2);
	cout<<"3: "<<VdotV(A,B)<<endl;
	if(VdotV(A,B)<=0.0){
		return false;
	}
	return true;

//	// The member result.sqrDistance is set in each block of the
//	// nested if-then-else statements. The remaining members are all
//	// set at the end of the function.
//
//	float const zero = static_cast<float>(0);
//	float const one = static_cast<float>(1);
//	float const two = static_cast<float>(2);
//
//
//	VmV(diff, trianglev0, point);
//
//
//	float a00 = VdotV(edge0, edge0);
//	float a01 = VdotV(edge0, edge1);
//	float a11 = VdotV(edge1, edge1);
//	float b0 = VdotV(diff, edge0);
//	float b1 = VdotV(diff, edge1);
//	float det = std::max(a00 * a11 - a01 * a01, (float)0.0);
//	float s = a01 * b1 - a11 * b0;
//	float t = a01 * b0 - a00 * b1;
//
//	// some bad triangles
//	if(a00==0.0||a01==0.0||a11==0.0||det==0.0){
//		return false;
//	}
//
//	if(a00==0.0||a01==0.0||a11==0.0||det==0.0){
//		printf("%f %f %f, %f %f %f, %f %f %f\n",
//				*(triangle+0)
//				,*(triangle+1)
//				,*(triangle+2)
//				,*(triangle+3)
//				,*(triangle+4)
//				,*(triangle+5)
//				,*(triangle+6)
//				,*(triangle+7)
//				,*(triangle+8));
//		printf("%f %f %f %f\n", a00, a01, a11, a00 * a11 - a01 * a01);
//	}
//
//	return s + t <= det && s>=zero && t>=zero;
}


}
