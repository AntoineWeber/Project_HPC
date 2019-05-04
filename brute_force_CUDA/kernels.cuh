
#ifndef __KERNELS_H__
#define __KERNELS_H__

__global__ void initialize_particles_uni(float *x_pos, float *y_pos, float *x_vel, float *y_vel, float *x_acc, float *y_acc, float *mass);
__global__ void initialize_particles_circle(float *x_pos, float *y_pos, float *x_vel, float *y_vel, float *x_acc, float *y_acc, float *mass);
__global__ void compute_forces(float *x_pos, float *y_pos, float *x_vel, float *y_vel, float *x_acc, float *y_acc, float *mass);

#endif