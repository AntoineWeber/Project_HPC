
#ifndef __KERNELS_H__
#define __KERNELS_H__

#define THETA 0.1
#define N_ITERATIONS 50
#define N_PARTICULES 1024
#define CHILD 4
#define G 6.67408e-11
#define EPSILON 0.001
#define MASS 100.0

#define BOUNDS 5
#define TIMESTEP 100
#define CIRCLE_OFFSET 3
#define FAR_SPACE 20


__global__ void initialize_particles_uni(double *x_pos, double *y_pos, double *x_vel, double *y_vel, double *x_acc, double *y_acc, double *mass);
__global__ void initialize_particles_circle(double *x_pos, double *y_pos, double *x_vel, double *y_vel, double *x_acc, double *y_acc, double *mass);

#endif