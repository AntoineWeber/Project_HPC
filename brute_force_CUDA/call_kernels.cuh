

#ifndef __CALL_KERNELS_CUH__
#define __CALL_KERNELS_CUH__
#include <cuda.h>
#include <cuda_runtime_api.h>


void initializeParticles(float *x_pos, float *y_pos, float *x_vel, float *y_vel, float *x_acc, float *y_acc,float *mass, dim3 gridSize, dim3 blockSize);
void computeForces(float *x_pos, float *y_pos, float *x_vel, float *y_vel, float *x_acc, float *y_acc, float *mass, dim3 gridSize, dim3 blockSize);


#endif