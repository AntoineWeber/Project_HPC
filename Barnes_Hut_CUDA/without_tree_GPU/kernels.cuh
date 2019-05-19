
#ifndef __KERNELS_H__
#define __KERNELS_H__
#include <stdio.h>
#include <math.h>
#include <curand.h>
#include <curand_kernel.h>
#include <cuda.h>

__global__ void initialize_particles_uni(double *x_pos, double *y_pos, double *x_vel, double *y_vel, double *x_acc, double *y_acc, double *mass);
__global__ void initialize_particles_circle(double *x_pos, double *y_pos, double *x_vel, double *y_vel, double *x_acc, double *y_acc, double *mass);

#endif