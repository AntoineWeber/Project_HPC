

#ifndef __CALL_KERNELS_CUH__
#define __CALL_KERNELS_CUH__
#include "kernels.cuh"
#include "barnesHut.hpp"
#include <cuda.h>
#include <cuda_runtime_api.h>


void initializeParticlesUni(double *x_pos, double *y_pos, double *x_vel, double *y_vel, double *x_acc, double *y_acc,double *mass, dim3 gridSize, dim3 blockSize);
void initializeParticlesCircle(double *x_pos, double *y_pos, double *x_vel, double *y_vel, double *x_acc, double *y_acc,double *mass, dim3 gridSize, dim3 blockSize);
void initializeParticles2Circles(double *x_pos, double *y_pos, double *x_vel, double *y_vel, double *x_acc, double *y_acc, double *mass, dim3 gridSize, dim3 blockSize);
void computeDisplacements(Node *d_tree, double* d_x, double* d_y, double* d_vx, double* d_vy, double* d_ax, double* d_ay, double* d_mass, dim3 gridSize, dim3 blockSize);


#endif