
#ifndef __KERNELS_H__
#define __KERNELS_H__
#include <stdio.h>
#include <math.h>
#include <curand.h>
#include <curand_kernel.h>
#include <cuda.h>
#include "barnesHut.hpp"

__global__ void initialize_particles_uni(double *x_pos, double *y_pos, double *x_vel, double *y_vel, double *x_acc, double *y_acc, double *mass);
__global__ void initialize_particles_circle(double *x_pos, double *y_pos, double *x_vel, double *y_vel, double *x_acc, double *y_acc, double *mass);
__global__ void compute_displacements(Node *d_tree, double* d_x, double* d_y, double* d_vx, double* d_vy, double* d_ax, double* d_ay, double* d_mass);
__global__ void compute_branch(int absOff, int depthOff, int nnode, Node *tree, double x, double y, double m, double &fx, double &fy);
#endif