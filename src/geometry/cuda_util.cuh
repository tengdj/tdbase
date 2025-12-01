/*
 * cuda_util.h
 *
 *  Created on: Dec 9, 2019
 *      Author: teng
 *
 *  can only be included by .cu files
 */

#ifndef MYCUDA_H_
#define MYCUDA_H_
#include "util.h"
#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>

namespace tdbase {

#define CUDA_SAFE_CALL(call)                                                                                                                         \
    do {                                                                                                                                             \
        cudaError_t err = call;                                                                                                                      \
        if (cudaSuccess != err) {                                                                                                                    \
            fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n", __FILE__, __LINE__, cudaGetErrorString(err));                              \
            exit(EXIT_FAILURE);                                                                                                                      \
        }                                                                                                                                            \
    } while (0);

inline void check_execution() {
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        log(cudaGetErrorString(err));
    }
}

// copy
__device__ inline void VcV_d(float Vr[3], const float V[3]) {
    Vr[0] = V[0];
    Vr[1] = V[1];
    Vr[2] = V[2];
}

// minus
__device__ inline void VmV_d(float Vr[3], const float V1[3], const float V2[3]) {
    Vr[0] = V1[0] - V2[0];
    Vr[1] = V1[1] - V2[1];
    Vr[2] = V1[2] - V2[2];
}

// plus
__device__ inline void VpV_d(float Vr[3], const float V1[3], const float V2[3]) {
    Vr[0] = V1[0] + V2[0];
    Vr[1] = V1[1] + V2[1];
    Vr[2] = V1[2] + V2[2];
}

// plus after product
__device__ inline void VpVxS_d(float Vr[3], const float V1[3], const float V2[3], float s) {
    Vr[0] = V1[0] + V2[0] * s;
    Vr[1] = V1[1] + V2[1] * s;
    Vr[2] = V1[2] + V2[2] * s;
}

// dot product
__device__ inline float VdotV_d(const float V1[3], const float V2[3]) {
    return (V1[0] * V2[0] + V1[1] * V2[1] + V1[2] * V2[2]);
}

__device__ inline void VcrossV_d(float Vr[3], const float V1[3], const float V2[3]) {
    Vr[0] = V1[1] * V2[2] - V1[2] * V2[1];
    Vr[1] = V1[2] * V2[0] - V1[0] * V2[2];
    Vr[2] = V1[0] * V2[1] - V1[1] * V2[0];
}

__device__ inline void sVpsV_2_d(float* Vr, float s1, const float* V1, float s2, const float* V2) {
    Vr[0] = s1 * V1[0] + s2 * V2[0];
    Vr[1] = s1 * V1[1] + s2 * V2[1];
}

__device__ inline void VxS_d(float Vr[3], const float V[3], float s) {
    Vr[0] = V[0] * s;
    Vr[1] = V[1] * s;
    Vr[2] = V[2] * s;
}

// Euclid distance
__device__ inline float VdistV2_d(const float V1[3], const float V2[3]) {
    return ((V1[0] - V2[0]) * (V1[0] - V2[0]) + (V1[1] - V2[1]) * (V1[1] - V2[1]) + (V1[2] - V2[2]) * (V1[2] - V2[2]));
}

}  // namespace tdbase

#endif /* MYCUDA_H_ */
