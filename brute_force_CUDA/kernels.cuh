
#ifndef __KERNELS_H__
#define __KERNELS_H__

#define PARTICLE_MASS 1
#define GRID_SIZE 32
#define BLOCK_SIZE 32

#define GRID_MIN -5
#define GRID_MAX 5

__global__ void initialize_particles(float *x_pos, float *y_pos, float *x_vel, float *y_vel, float *x_acc, float *y_acc, float *mass, int nparticles);

#endif