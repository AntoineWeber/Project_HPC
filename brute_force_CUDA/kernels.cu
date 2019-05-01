
#include <stdio.h>
#include <stdlib.h>
#include "kernels.cuh"


__global__ void initialize_particles(float *x_pos, float *y_pos, float *x_vel, float *y_vel, float *x_acc, float *y_acc, float *mass, int nparticles)
{
    int max_rand = 32767;
    int i = threadIdx.x + blockIdx.x*blockDim.x;

    mass[i] = (float)PARTICLE_MASS;
    x[i] = (rand() % (GRID_MAX - GRID_MIN + 1)) + GRID_MIN;
    y[i] = (rand() % (GRID_MAX - GRID_MIN + 1)) + GRID_MIN;

    // set velocity to 0
    x_vel[i] = 0;
    y_vel[i] = 0;

    // set acceleration to zero
    x_acc[i] = 0.0;
    y_acc[i] = 0.0;
}