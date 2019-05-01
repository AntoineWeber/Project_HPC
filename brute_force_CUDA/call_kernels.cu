


#include "call_kernels.cuh"
#include "kernels.cuh"
#include "brute_force.h"


void initializeParticles(float *x_pos, float *y_pos, float *x_vel, float *y_vel, float *x_acc, float *y_acc, float *d_mass)
{
    dim3 gridSize(GRID_SIZE);
    dim3 blockSize(BLOCK_SIZE);
	initialize_particles<<< gridSize, blockSize>>>(x_pos, y_pos, x_vel, y_vel, x_acc, y_acc, d_mass, N_PARTICLES);//, d_state);
}

