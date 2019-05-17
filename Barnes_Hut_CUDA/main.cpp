


#include <iostream>
#include <fstream>
#include <vector>
#include <iostream>
#include <memory>
#include <random>
#include <limits>
#include <cmath>

#include "kernels.cuh"

#include <cuda.h>
#include <cuda_runtime_api.h>

//#define SAVE true


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
    // initialize time for TTS (time to solution)
    Timer t1;

    // define grid and block size
    dim3 gridSize(GRID_SIZE);
    dim3 blockSize(BLOCK_SIZE);

    initiateDevice();

    // allocate memory for the particle array
    unsigned int size_particles = N_PARTICLES; // plane coordinates
    unsigned int mem_size_particles = sizeof(double) * size_particles;

    double* d_pos_x;
    double* d_pos_y;
    double* d_vel_x;
    double* d_vel_y;
    double* d_acc_x;
    double* d_acc_y;
    double* d_mass;

    gpuErrchk(cudaMalloc((void**) &d_pos_x, mem_size_particles));
    gpuErrchk(cudaMalloc((void**) &d_pos_y, mem_size_particles));
    gpuErrchk(cudaMalloc((void**) &d_vel_x, mem_size_particles));
    gpuErrchk(cudaMalloc((void**) &d_vel_y, mem_size_particles));
    gpuErrchk(cudaMalloc((void**) &d_acc_x, mem_size_particles));
    gpuErrchk(cudaMalloc((void**) &d_acc_y, mem_size_particles));
    gpuErrchk(cudaMalloc((void**) &d_mass, mem_size_particles));

    //initializeParticlesUni(d_pos_x, d_pos_y, d_vel_x, d_vel_y, d_acc_x, d_acc_y, d_mass, gridSize, blockSize);
    initializeParticlesCircle(d_pos_x, d_pos_y, d_vel_x, d_vel_y, d_acc_x, d_acc_y, d_mass, gridSize, blockSize);

    // used for output saving
    #ifdef SAVE
        double *h_pos_x = (double*)malloc(mem_size_particles);
        double *h_pos_y = (double*)malloc(mem_size_particles);
        std::vector<Position> output(N_PARTICLES);

        //open file
        std::ofstream myFile("trajectories.txt");

        cudaMemcpy(h_pos_x, d_pos_x, mem_size_particles, cudaMemcpyDeviceToHost);
        cudaMemcpy(h_pos_y, d_pos_y, mem_size_particles, cudaMemcpyDeviceToHost);

        output = fillArray(h_pos_x, h_pos_y);
        for(std::vector<Position>::const_iterator i = output.begin(); i != output.end(); ++i)
        {
            myFile << (*i).x << " " << (*i).y << "\n";
        }
    #endif

    double elapsed_ini = t1.elapsed();
    t1.reset();

    for (unsigned int iter=0; iter<ITERATIONS; iter++)
    {
        /*
        resetTree()
        buildTree()
        computeDisplacement()
        */
        #ifdef SAVE
            cudaMemcpy(h_pos_x, d_pos_x, mem_size_particles, cudaMemcpyDeviceToHost);
            cudaMemcpy(h_pos_y, d_pos_y, mem_size_particles, cudaMemcpyDeviceToHost);

            output = fillArray(h_pos_x, h_pos_y);
            for(std::vector<Position>::const_iterator i = output.begin(); i != output.end(); ++i)
            {
                myFile << (*i).x << " " << (*i).y << "\n";
            }
        #endif
    }
    double elapsed_compute = t1.elapsed();


    // GTX 1070 has a theoretical 6.5 TFLOPS
    std::cout << std::endl;
    std::cout << "Elapsed time for initialization : " << elapsed_ini << " s" << std::endl;
    std::cout << "Elapsed time for computation    : " << elapsed_compute << " s" << std::endl;
    std::cout << "computing " << ITERATIONS << " steps with " << N_PARTICLES << " particles." <<std::endl;
    std::cout << std::endl;

    gpuErrchk(cudaFree(d_pos_x));
    gpuErrchk(cudaFree(d_pos_y));
    gpuErrchk(cudaFree(d_vel_x));
    gpuErrchk(cudaFree(d_vel_y));
    gpuErrchk(cudaFree(d_acc_x));
    gpuErrchk(cudaFree(d_acc_y));
    gpuErrchk(cudaFree(d_mass));

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

#ifdef SAVE
std::vector<Position> fillArray(double *pos_x, double *pos_y)
{
    std::vector<Position> output(N_PARTICLES);
    for (int i=0; i<N_PARTICLES; i++)
    {
        output[i].x = pos_x[i];
        output[i].y = pos_y[i];
    }
    return output;
}
#endif
