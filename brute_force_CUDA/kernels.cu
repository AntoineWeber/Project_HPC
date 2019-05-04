
#include <stdio.h>
#include <math.h>
#include <curand.h>
#include <curand_kernel.h>
#include <cuda.h>
#include "kernels.cuh"
#include "brute_force.h"


__global__ void initialize_particles_uni(float *x_pos, float *y_pos, float *x_vel, float *y_vel, float *x_acc, float *y_acc, float *mass)
{
    int i = threadIdx.x + blockIdx.x*blockDim.x;
    int stride = gridDim.x * blockDim.x;
    int offset = 0;

    // to initialize the cuda rand
    curandState state;
    curand_init(clock64(), i, 0, &state);

    while (i + offset < N_PARTICLES)
    {
        mass[i] = (float)PARTICLE_MASS;

        x_pos[i] = curand_uniform(&state)*GRID_MAX*2 + GRID_MIN;
        y_pos[i] = curand_uniform(&state)*GRID_MAX*2 + GRID_MIN;

        // set velocity to 0
        x_vel[i] = 0.0;
        y_vel[i] = 0.0;

        // set acceleration to 0
        x_acc[i] = 0.0;
        y_acc[i] = 0.0;

        offset += stride;
    }
}

__global__ void initialize_particles_circle(float *x_pos, float *y_pos, float *x_vel, float *y_vel, float *x_acc, float *y_acc, float *mass)
{
    int i = threadIdx.x + blockIdx.x*blockDim.x;
    int stride = gridDim.x * blockDim.x;
    int offset = 0;

    // to initialize the cuda rand
    curandState state;
    curand_init(clock64(), i, 0, &state);

    while (i + offset < N_PARTICLES)
    {
        mass[i] = (float)PARTICLE_MASS;

        float r = R_OFFSET + curand_uniform(&state)*R_MAX;
        float alpha = curand_uniform(&state)*2*PI;

        x_pos[i] = r*cos(alpha);
        y_pos[i] = r*sin(alpha);

        // set velocity to 0
        x_vel[i] = 0.0;
        y_vel[i] = 0.0;

        // set acceleration to 0
        x_acc[i] = 0.0;
        y_acc[i] = 0.0;

        offset += stride;
    }

}



__global__ void compute_forces(float *x_pos, float *y_pos, float *x_vel, float *y_vel, float *x_acc, float *y_acc, float *mass)
{
    int i = threadIdx.x + blockIdx.x*blockDim.x;
    int stride = gridDim.x * blockDim.x;
    int offset = 0;

    float fx=0;
    float fy=0;
    // loop on given particles if not enough threads
    while (i + offset < N_PARTICLES)
    {
        #pragma unroll
        for (int j=0; j<N_PARTICLES; j++)
        {
            if (j != i)
            {
                float r = sqrt((x_pos[j]-x_pos[i+offset])*(x_pos[j]-x_pos[i+offset]) + (y_pos[j]-y_pos[i+offset])*(y_pos[j]-y_pos[i+offset]));
                // may result in a division by 0 otherwise
                // should do something at that point but not sure how to handle it
                if (r > 0.001)
                {
                    fx += (G*mass[i+offset]*mass[j]*(x_pos[j]-x_pos[i+offset]))/(sqrt((x_pos[j]-x_pos[i+offset])*(x_pos[j]-x_pos[i+offset])));
                    fy += (G*mass[i+offset]*mass[j]*(y_pos[j]-y_pos[i+offset]))/(sqrt((y_pos[j]-y_pos[i+offset])*(y_pos[j]-y_pos[i+offset])));
                }
            }
        }

        // F = ma -> a = F/m
        x_acc[i+offset] = fx / mass[i+offset];
        y_acc[i+offset] = fy / mass[i+offset];
        
        offset += stride;
    }
    __syncthreads();
    
    // have to define this second loop otherwise would move particles before having computed the force for all of them
    offset = 0;
    while (i + offset < N_PARTICLES)
    {
        x_pos[i+offset] += 0.5*x_acc[i+offset]*TIMESTEP*TIMESTEP + x_vel[i+offset]*TIMESTEP;
        y_pos[i+offset] += 0.5*y_acc[i+offset]*TIMESTEP*TIMESTEP + y_vel[i+offset]*TIMESTEP;

        x_vel[i+offset] += x_acc[i+offset]*TIMESTEP;
        y_vel[i+offset] += y_acc[i+offset]*TIMESTEP;

        offset+=stride;
    }
}
