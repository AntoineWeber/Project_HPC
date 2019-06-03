#include "call_kernels.cuh"

// File used to call CUDA kernels from C++ code

void initializeParticlesUni(double *x_pos, double *y_pos, double *x_vel, double *y_vel, double *x_acc, double *y_acc, double *mass, dim3 gridSize, dim3 blockSize)
{
	initialize_particles_uni<<< gridSize, blockSize>>>(x_pos, y_pos, x_vel, y_vel, x_acc, y_acc, mass);
}

void initializeParticlesCircle(double *x_pos, double *y_pos, double *x_vel, double *y_vel, double *x_acc, double *y_acc, double *mass, dim3 gridSize, dim3 blockSize)
{
	initialize_particles_circle<<< gridSize, blockSize>>>(x_pos, y_pos, x_vel, y_vel, x_acc, y_acc, mass);
}
void initializeParticles2Circles(double *x_pos, double *y_pos, double *x_vel, double *y_vel, double *x_acc, double *y_acc, double *mass, dim3 gridSize, dim3 blockSize)
{
	initialize_particles_2_circles<<< gridSize, blockSize>>>(x_pos, y_pos, x_vel, y_vel, x_acc, y_acc, mass);
}

void computeDisplacements(Node *d_tree, double* d_x, double* d_y, double* d_vx, double* d_vy, double* d_ax, double* d_ay, double* d_mass, dim3 gridSize, dim3 blockSize)
{
	compute_displacements<<<gridSize,blockSize>>>(d_tree, d_x, d_y, d_vx, d_vy, d_ax, d_ay, d_mass);
}
