
#include <stdio.h>
#include <math.h>
#include <curand.h>
#include <curand_kernel.h>
#include <cuda.h>
#include "kernels.cuh"
#include "brute_force.h"


__global__ void initialize_particles_uni(double *x_pos, double *y_pos, double *x_vel, double *y_vel, double *x_acc, double *y_acc, double *mass)
{
    int i = threadIdx.x + blockIdx.x*blockDim.x;
    int stride = gridDim.x * blockDim.x;
    int offset = 0;

    // to initialize the cuda rand
    curandState state;
    curand_init(clock64(), i, 0, &state);

    while (i + offset < N_PARTICLES)
    {
        mass[i + offset] = (double)PARTICLE_MASS;

        x_pos[i + offset] = curand_uniform(&state)*GRID_MAX*2 + GRID_MIN;
        y_pos[i + offset] = curand_uniform(&state)*GRID_MAX*2 + GRID_MIN;

        // set velocity to 0
        x_vel[i + offset] = 0.0;
        y_vel[i + offset] = 0.0;

        // set acceleration to 0
        x_acc[i + offset] = 0.0;
        y_acc[i + offset] = 0.0;

        offset += stride;
    }
}

__global__ void initialize_particles_circle(double *x_pos, double *y_pos, double *x_vel, double *y_vel, double *x_acc, double *y_acc, double *mass)
{
    int i = threadIdx.x + blockIdx.x*blockDim.x;
    int stride = gridDim.x * blockDim.x;
    int offset = 0;

    // to initialize the cuda rand
    curandState state;
    curand_init(clock64(), i, 0, &state);

    while (i + offset < N_PARTICLES)
    {
        mass[i + offset] = (double)PARTICLE_MASS;

        double r = R_OFFSET + curand_uniform(&state)*R_MAX;
        double alpha = curand_uniform(&state)*2*M_PI;

        x_pos[i + offset] = r*cos(alpha);
        y_pos[i + offset] = r*sin(alpha);

        // set velocity to 0
        x_vel[i + offset] = 0.0;
        y_vel[i + offset] = 0.0;

        // set acceleration to 0
        x_acc[i + offset] = 0.0;
        y_acc[i + offset] = 0.0;

        offset += stride;
    }

}
