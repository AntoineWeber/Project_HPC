
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
void copyTreetoCuda(QuadTree* node, d_Node* d_tree);

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

    dim3 gridSize(GRID_SIZE);
    dim3 blockSize(BLOCK_SIZE);

    // used for output saving
    #ifdef SAVE
        //open file
        std::ofstream myFile("trajectories.txt");
    #endif

    t1.reset();
    Particles allParticles;

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

    // allocate quadtree for device
    d_Node to_copy;
    d_Node *d_tree;
    gpuErrchk(cudaMalloc((void **)&d_tree, sizeof(d_Node)));
    gpuErrchk(cudaMemcpy(d_tree, &to_copy, sizeof(d_Node), cudaMemcpyHostToDevice));

    double elapsed_ini = t1.elapsed();

    std::cout << "NOT WORKING. COPY TOO SLOW AND MOREOVER NOT CORRECT. ABORTING" << std::endl;
    
    #ifdef SAVE
        allParticles.saveToFile(&myFile);
    #endif

    t1.reset();
    for (unsigned int i=0; i<N_ITERATIONS; i++)
    {
        std::cout << "iteration " << i << std::endl;
        allParticles.resetTree();
        allParticles.buildTree();

        std::cout << "entering copy on GPU" << std::endl;
        copyTreetoCuda(&(allParticles.m_tree), d_tree);
        std::cout << "leaving copy on GPU" << std::endl;

        #ifdef SAVE
            allParticles.saveToFile(&myFile);
        #endif
    }
    double elapsed_compute = t1.elapsed();
    allParticles.resetTree();

    std::cout << "The initialisation took : " << elapsed_ini << "seconds." << std::endl;
    std::cout << "While the computation of the " << N_ITERATIONS << " steps took " << elapsed_compute << "seconds." << std::endl;
    std::cout << "with " << N_PARTICULES << " particles." << std::endl;

    free(d_tree);

    return 1;
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

void copyTreetoCuda(QuadTree* node, d_Node* d_tree)
{
    d_Node current;
    // not working as i want it to be which is normal
    d_tree = new d_Node;

    for (unsigned int i=0; i<node->m_children.size(); i++)
    {
        if (node->m_children[i] != nullptr)
        {
            // allocate device node
            gpuErrchk(cudaMalloc((void **) &(d_tree->m_children[i]), sizeof(d_Node)));
            // copy value in device class
            current.hasChildren = node->m_children[i]->hasChildren;
            current.m_av_mass = node->m_children[i]->m_av_mass;
            current.m_s = node->m_children[i]->m_s;
            current.m_x_center = node->m_children[i]->m_x_center;
            current.m_y_center = node->m_children[i]->m_y_center;

            // copy to GPU
            gpuErrchk(cudaMemcpy(d_tree->m_children[i], &current, sizeof(d_Node), cudaMemcpyHostToDevice));
            // recursive call
            if (node->m_children[i]->hasChildren)
            {
                copyTreetoCuda(node->m_children[i], d_tree->m_children[i]);   
            }
        }
    }
}


