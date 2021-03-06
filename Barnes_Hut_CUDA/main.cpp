
#include <chrono>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "barnesHut.hpp"
#include "call_kernels.cuh"

#include <cuda.h>
#include <cuda_runtime_api.h>
#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>

void initiateDevice();

// Timer class
class Timer
{
public:
    Timer() : beg_(clock_::now()) {}
    void reset() { beg_ = clock_::now(); }
    double elapsed() const {
        return std::chrono::duration_cast<second_>
            (clock_::now() - beg_).count(); }

private:
    typedef std::chrono::high_resolution_clock clock_;
    typedef std::chrono::duration<double, std::ratio<1> > second_;
    std::chrono::time_point<clock_> beg_;
};

// CUDA error checking
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}


////////////////////////////////////////////////////////////////////////////////
// Program main : 2D
////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
    Timer t1;

    initiateDevice();

    // set grid and block size
    dim3 gridSize(GRID_SIZE);
    dim3 blockSize(BLOCK_SIZE);

    // used for output saving
    #ifdef SAVE
        //open file
        std::ofstream myFile("trajectories.txt");
    #endif

    t1.reset();

    // compute the number of nodes in the tree considering 
    // MAX_DEPTH layers (starts at 0)
    int num_nodes = 0;
    for (unsigned int i = 0; i<MAX_DEPTH; i++)
    {
        num_nodes += pow(4.0,i);
    }
    std::cout << "generating " << num_nodes << " nodes." << std::endl;

    // needs num_nodes to allocate the host fixed-size tree
    Particles allParticles(num_nodes);

    // allocate tree for device
    Node *d_tree;
    gpuErrchk(cudaMalloc((void**) &d_tree, num_nodes*sizeof(Node)));

    // copy vectors to GPU
    thrust::device_vector<double> d_x  = allParticles.m_x;
    thrust::device_vector<double> d_y  = allParticles.m_y;
    thrust::device_vector<double> d_vx = allParticles.m_vx;
    thrust::device_vector<double> d_vy = allParticles.m_vy;
    thrust::device_vector<double> d_ax = allParticles.m_ax;
    thrust::device_vector<double> d_ay = allParticles.m_ay;
    thrust::device_vector<double> d_mass = allParticles.m_mass;

    // initialize positions of particles
    initializeParticlesCircle(thrust::raw_pointer_cast(&d_x[0]), thrust::raw_pointer_cast(&d_y[0]), 
                              thrust::raw_pointer_cast(&d_vx[0]), thrust::raw_pointer_cast(&d_vy[0]),
                              thrust::raw_pointer_cast(&d_ax[0]), thrust::raw_pointer_cast(&d_ay[0]), 
                              thrust::raw_pointer_cast(&d_mass[0]), gridSize, blockSize);

    // copy position to host
    allParticles.m_x = d_x;
    allParticles.m_y = d_y;
    allParticles.m_vx = d_vx;
    allParticles.m_vy = d_vy;
    allParticles.m_ax = d_ax;
    allParticles.m_ay = d_ay;
    allParticles.m_mass = d_mass;

    double elapsed_ini = t1.elapsed();
    
    #ifdef SAVE
        allParticles.saveToFile(&myFile);
    #endif

    t1.reset();
    for (unsigned int i=0; i<N_ITERATIONS; i++)
    {
        // Clean the tree and allocate another one
        allParticles.resetTree(num_nodes);
        // Build the tree locally
        allParticles.buildTree();
        // Copy the tree to device
        gpuErrchk(cudaMemcpy(d_tree, &allParticles.m_tree[0], num_nodes*sizeof(Node), cudaMemcpyHostToDevice));
        // Compute displacement
        computeDisplacements(d_tree, thrust::raw_pointer_cast(&d_x[0]),
                             thrust::raw_pointer_cast(&d_y[0]), thrust::raw_pointer_cast(&d_vx[0]),
                             thrust::raw_pointer_cast(&d_vy[0]), thrust::raw_pointer_cast(&d_ax[0]),
                             thrust::raw_pointer_cast(&d_ay[0]), thrust::raw_pointer_cast(&d_mass[0]), gridSize, blockSize);
        // synchronize before new loop
        cudaDeviceSynchronize();
        
        // copy particles position back to local to update tree
        allParticles.m_x = d_x;
        allParticles.m_y = d_y;
                                                
        #ifdef SAVE
            allParticles.saveToFile(&myFile);
        #endif
    }

    // free pointer for device tree
    if (d_tree != nullptr)
    {
        gpuErrchk(cudaFree(d_tree));
        d_tree = nullptr;
    }

    double elapsed_compute = t1.elapsed();

    std::cout << "The initialisation took : " << elapsed_ini << "seconds." << std::endl;
    std::cout << "While the computation of the " << N_ITERATIONS << " steps took " << elapsed_compute << "seconds." << std::endl;
    std::cout << "with " << N_PARTICULES << " particles." << std::endl;

    return 0;
}


// Function checking if the GPU can be used for computation
void initiateDevice()
{
    // By default, we use device 0, otherwise we override the device ID based on what is provided at the command line
    int devID = 0;

    cudaError_t error;
    cudaDeviceProp deviceProp;
    error = cudaGetDevice(&devID);


    error = cudaGetDeviceProperties(&deviceProp, devID);
    std::cout << "Memory bus width : " << deviceProp.memoryBusWidth << std::endl;
    std::cout << "Memory clock rate : " << deviceProp.memoryClockRate << std::endl;

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
