
// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <math.h>

#include "brute_force.h"
#include "call_kernels.cuh"
#include <cuda.h>
#include <cuda_runtime_api.h>


////////////////////////////////////////////////////////////////////////////////
// Program main : 2D
////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
    initiateDevice();

    // allocate memory for the particle array
    unsigned int size_particles = N_PARTICLES * 2; // plane coordinates
    unsigned int mem_size_particles = sizeof(float) * size_particles;

    float* d_pos_x;
    float* d_pos_y;
    float* d_vel_x;
    float* d_vel_y;
    float* d_acc_x;
    float* d_acc_y;
    float* d_mass;

    dim3 gridSize(GRID_SIZE);
    dim3 blockSize(BLOCK_SIZE);

    cudaMalloc((void**) &d_pos_x, mem_size_particles);
    cudaMalloc((void**) &d_pos_y, mem_size_particles);
    cudaMalloc((void**) &d_vel_x, mem_size_particles);
    cudaMalloc((void**) &d_vel_y, mem_size_particles);
    cudaMalloc((void**) &d_acc_x, mem_size_particles);
    cudaMalloc((void**) &d_acc_y, mem_size_particles);
    cudaMalloc((void**) &d_mass, mem_size_particles);


    Timer t1;
    initializeParticles(d_pos_x, d_pos_y, d_vel_x, d_vel_y, d_acc_x, d_acc_y, d_mass, gridSize, blockSize);

    /*
    float *h_pos_x = (float*)malloc(mem_size_particles);
    float *post_h_pos_x = (float*)malloc(mem_size_particles);
    cudaMemcpy(h_pos_x, d_pos_x, mem_size_particles, cudaMemcpyDeviceToHost);
    */

    for (unsigned int iter=0; iter<ITERATIONS; iter++)
    {
        computeForces(d_pos_x, d_pos_y, d_vel_x, d_vel_y, d_acc_x, d_acc_y, d_mass, gridSize, blockSize);
        //cudaMemcpy(post_h_pos_x, d_pos_x, mem_size_particles, cudaMemcpyDeviceToHost);
    }
    double elapsed = t1.elapsed();


    // GTX 1070 has a theoretical 6.5 TFLOPS
    std::cout << std::endl;
    std::cout << "Elapsed time : " << elapsed << " s" << std::endl;
    std::cout << "computing " << ITERATIONS << " steps with " << N_PARTICLES << " particles." <<std::endl;
    std::cout << std::endl;

    /*
    for (unsigned int j=0; j<N_PARTICLES; j++)
    {
        std::cout << post_h_pos_x[j] - h_pos_x[j] << " ";
    }
    std::cout << std::endl;
    */

    cudaFree(d_pos_x);
    cudaFree(d_pos_y);
    cudaFree(d_vel_x);
    cudaFree(d_vel_y);
    cudaFree(d_acc_x);
    cudaFree(d_acc_y);
    cudaFree(d_mass);

    exit(EXIT_SUCCESS);
}

void initiateDevice()
{
    // By default, we use device 0, otherwise we override the device ID based on what is provided at the command line
    int devID = 0;

    cudaError_t error;
    cudaDeviceProp deviceProp;
    error = cudaGetDevice(&devID);


    error = cudaGetDeviceProperties(&deviceProp, devID);

    if (deviceProp.computeMode == cudaComputeModeProhibited)
    {
        fprintf(stderr, "Error: device is running in <Compute Mode Prohibited>, no threads can use ::cudaSetDevice().\n");
        exit(EXIT_SUCCESS);
    }

    if (error != cudaSuccess)
    {
        printf("cudaGetDeviceProperties returned error code %d, line(%d)\n", error, __LINE__);
    }
    else
    {
        printf("GPU Device %d: \"%s\" with compute capability %d.%d\n\n", devID, deviceProp.name, deviceProp.major, deviceProp.minor);
    }
}