#ifndef HISPEED_GEOMETRY_H
#define HISPEED_GEOMETRY_H

#include <math.h>
#include <stdio.h>
#include <iostream>
#include <float.h>

#include "../util/util.h"
using namespace std;

namespace hispeed{

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

struct dist_param{
	int id;
	const float *data1;
	const float *data2;
	const float *result;
	int size1;
	int size2;
	float dist;
};

float TriDist(const float *S, const float *T);
float TriDist_batch(const float *S, const float * T, size_t s1, size_t s2, int num_threads);

float SegDist_single(const float *data1, const float *data2, size_t size1, size_t size2);
float SegDist_batch(const float *S1, const float *S2, size_t size1, size_t size2, int num_threads);
void SegDist_batch_gpu(const float *S, const float *T, int batch_size, int batch_num, float *result);
}
#endif
