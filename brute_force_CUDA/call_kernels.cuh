

#ifndef __CALL_KERNELS_CUH__
#define __CALL_KERNELS_CUH__
#include <cuda.h>
#include <cuda_runtime_api.h>


void initializeParticlesUni(double *x_pos, double *y_pos, double *x_vel, double *y_vel, double *x_acc, double *y_acc,double *mass, dim3 gridSize, dim3 blockSize);
void initializeParticlesCircle(double *x_pos, double *y_pos, double *x_vel, double *y_vel, double *x_acc, double *y_acc,double *mass, dim3 gridSize, dim3 blockSize);
void computeForces(double *x_pos, double *y_pos, double *x_vel, double *y_vel, double *x_acc, double *y_acc, double *mass, dim3 gridSize, dim3 blockSize);


#endif