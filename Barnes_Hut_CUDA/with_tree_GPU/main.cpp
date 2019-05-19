#include <iostream>
#include <fstream>
#include <vector>
#include <iostream>
#include <memory>
#include <random>
#include <limits>
#include <cmath>
#include <chrono>

#include "call_kernels.cuh"
#include "kernels.cuh"

void initiateDevice();


// small timer class
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
    // initialize time for TTS (time to solution)
    Timer t1;

    initiateDevice();

    // define grid and block size
    dim3 gridSize(GRID_SIZE);
    dim3 blockSize(BLOCK_SIZE);

    // instanciate particles class on GPU
    Particles *allParticles;
    gpuErrchk(cudaMalloc((void **) &allParticles, N_PARTICLES*sizeof(Particles)));

    // instanciate the tree on GPU. Defined a maximum depth to have a fixed size.
    QuadTree *allNodes;

    // how many nodes in the tree
    int num_nodes = 0;
    for (unsigned int i = 0; i<MAX_DEPTH; i++)
    {
        num_nodes += pow(4.0,i);
    }
    gpuErrchk(cudaMalloc((void **) &allNodes, num_nodes * sizeof(QuadTree)))
    // now set first node to get all particles
    QuadTree first_node;
    first_node.m_begin_particle = 0;
    first_node.m_end_particle = N_PARTICLES - 1;
    cudaMemcpy(allNodes, &first_node, sizeof(QuadTree), cudaMemcpyHostToDevice);

    // initialize particle positions
    initializeParticles(allParticles, gridSize, blockSize);
    gpuErrchk(cudaDeviceSynchronize());

    // used for output saving
    #ifdef SAVE
        Particles *h_partic = new Particles[N_PARTICLES];

        //open file
        std::ofstream myFile("trajectories.txt");

        cudaMemcpy(h_partic, allParticles, N_PARTICLES*sizeof(Particles), cudaMemcpyDeviceToHost);

        for(unsigned int i = 0; i<N_PARTICLES; i++)
        {
            myFile << h_partic[i].m_x << " " << h_partic[i].m_y << "\n";
        }
    #endif

    double elapsed_ini = t1.elapsed();
    t1.reset();


    for (unsigned int iter=0; iter<N_ITERATIONS; iter++)
    {
        /*
        resetTree()
        buildTree()
        computeDisplacement()
        */
        #ifdef SAVE
            cudaMemcpy(h_partic, allParticles, N_PARTICLES*sizeof(Particles), cudaMemcpyDeviceToHost);

            for(unsigned int i = 0; i<N_PARTICLES; i++)
            {
                myFile << h_partic[i].m_x << " " << h_partic[i].m_y << "\n";
            }
        #endif
    }
    double elapsed_compute = t1.elapsed();


    // GTX 1070 has a theoretical 6.5 TFLOPS
    std::cout << std::endl;
    std::cout << "Elapsed time for initialization : " << elapsed_ini << " s" << std::endl;
    std::cout << "Elapsed time for computation    : " << elapsed_compute << " s" << std::endl;
    std::cout << "computing " << N_ITERATIONS << " steps with " << N_PARTICLES << " particles." <<std::endl;
    std::cout << std::endl;
    
    #ifdef SAVE
        free(h_partic);
    #endif


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
