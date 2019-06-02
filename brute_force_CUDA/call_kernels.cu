


#include "call_kernels.cuh"
#include "kernels.cuh"
#include "brute_force.h"

// File used to call the cuda kernels from the C++ file.

void initializeParticlesUni(double *x_pos, double *y_pos, double *x_vel, double *y_vel, double *x_acc, double *y_acc, double *mass, dim3 gridSize, dim3 blockSize)
{
	initialize_particles_uni<<< gridSize, blockSize>>>(x_pos, y_pos, x_vel, y_vel, x_acc, y_acc, mass);
}

void initializeParticlesCircle(double *x_pos, double *y_pos, double *x_vel, double *y_vel, double *x_acc, double *y_acc, double *mass, dim3 gridSize, dim3 blockSize)
{
	initialize_particles_circle<<< gridSize, blockSize>>>(x_pos, y_pos, x_vel, y_vel, x_acc, y_acc, mass);
}

void computeForces(double *x_pos, double *y_pos, double *x_vel, double *y_vel, double *x_acc, double *y_acc, double *mass, dim3 gridSize, dim3 blockSize)
{
    compute_forces<<<gridSize, blockSize>>>(x_pos, y_pos, x_vel, y_vel, x_acc, y_acc, mass);
}
