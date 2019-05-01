
#include <stdio.h>

#include <curand.h>
#include <curand_kernel.h>
#include <cuda.h>
#include "kernels.cuh"


__global__ void initialize_particles(float *x_pos, float *y_pos, float *x_vel, float *y_vel, float *x_acc, float *y_acc, float *mass, int nparticles)//, curandState *state)
{
    int i = threadIdx.x + blockIdx.x*blockDim.x;
    int stride = gridDim.x * blockDim.x;
    int offset = 0;

    // to initialize the cuda rand
    curandState state;
    curand_init(clock64(), i, 0, &state);

    while (i + offset < nparticles)
    {
        mass[i] = (float)PARTICLE_MASS;

        x_pos[i] = curand_uniform(&state)*GRID_MAX*2 + GRID_MIN;
        y_pos[i] = curand_uniform(&state)*GRID_MAX*2 + GRID_MIN;

        // set velocity to 0
        x_vel[i] = 0.0;
        y_vel[i] = 0.0;

        // set acceleration to zero
        x_acc[i] = 0.0;
        y_acc[i] = 0.0;

        offset += stride;
    }
}
