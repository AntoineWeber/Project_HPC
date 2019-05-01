
// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "brute_force.h"
#include "call_kernels.cuh"
#include <cuda.h>


////////////////////////////////////////////////////////////////////////////////
// Program main : 2D
////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
    initiateDevice();

    // allocate memory for the particle array
    unsigned int size_particles = N_PARTICLES * 2; // plane coordinates
    unsigned int mem_size_particles = sizeof(float) * size_particles;

    float* d_pos_x,d_pos_y,d_vel_x,d_vel_y,d_acc_x,d_acc_y,d_mass;
    cudaMalloc((void**) &d_pos_x, mem_size);
    cudaMalloc((void**) &d_pos_y, mem_size);
    cudaMalloc((void**) &d_vel_x, mem_size);
    cudaMalloc((void**) &d_vel_y, mem_size);
    cudaMalloc((void**) &d_acc_x, mem_size);
    cudaMalloc((void**) &d_acc_y, mem_size);
    cudaMalloc((void**) &d_mass, mem_size);

    initializeParticles(d_pos_x, d_pos_y, d_vel_x, d_vel_y, d_acc_x, d_acc_y, d_mass);

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