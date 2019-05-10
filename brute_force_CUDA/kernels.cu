
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



__global__ void compute_forces(double *x_pos, double *y_pos, double *x_vel, double *y_vel, double *x_acc, double *y_acc, double *mass)
{
    int i = threadIdx.x + blockIdx.x*blockDim.x;
    int stride = gridDim.x * blockDim.x;
    int offset = 0;

    double fx=0;
    double fy=0;
    // loop on given particles if not enough threads
    while (i + offset < N_PARTICLES)
    {
        fx = 0;
        fy = 0;
        // if outside the defined box, stop taking this particle into consideration
        if (x_pos[i + offset] > SIMU_BOUND_X || x_pos[i + offset] < -SIMU_BOUND_X || y_pos[i + offset] > SIMU_BOUND_Y || y_pos[i + offset] < -SIMU_BOUND_Y)
        {
            x_acc[i + offset] = 0;
            y_acc[i + offset] = 0;
        }
        else
        {
            // when looping on itself, will just add a 0 composant to the force.
            for (int j=0; j<N_PARTICLES; j++)
            {
                // if target in the defined grid
                if (!(x_pos[j] > SIMU_BOUND_X || x_pos[j] < -SIMU_BOUND_X || y_pos[j] > SIMU_BOUND_Y || y_pos[j] < -SIMU_BOUND_Y))
                {
                    double r = sqrt((x_pos[j]-x_pos[i+offset])*(x_pos[j]-x_pos[i+offset]) + (y_pos[j]-y_pos[i+offset])*(y_pos[j]-y_pos[i+offset]));
                    // if distance between the points smaller than epsilon, fix d = epsilon for computation.
                    if (r > EPSILON)
                    {
                        fx += (G*mass[i+offset]*mass[j]*(x_pos[j]-x_pos[i+offset]))/r;
                        fy += (G*mass[i+offset]*mass[j]*(y_pos[j]-y_pos[i+offset]))/r;
                    }
                    else
                    {
                        fx += (G*mass[i+offset]*mass[j]*(x_pos[j]-x_pos[i+offset]))/EPSILON;
                        fy += (G*mass[i+offset]*mass[j]*(y_pos[j]-y_pos[i+offset]))/EPSILON;
                    }
                }
            }

            // F = ma -> a = F/m
            x_acc[i+offset] = fx / mass[i+offset];
            y_acc[i+offset] = fy / mass[i+offset];
        }
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
