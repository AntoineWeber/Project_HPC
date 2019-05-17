

#ifndef __CALL_KERNELS_CUH__
#define __CALL_KERNELS_CUH__

#define THETA 0.3
#define N_ITERATIONS 100
#define CHILD 4
#define G 1.0
#define EPSILON 0.001
#define MASS 1.0

//#define SAVE true

#define BOUNDS 5
#define TIMESTEP 0.001
#define CIRCLE_OFFSET 3
#define FAR_SPACE 20

#define N_PARTICLES 2048
#define GRID_SIZE 32
#define BLOCK_SIZE 32
#define WARP_SIZE 32

// 3 is the margin that I take
#define MAX_DEPTH (int)(log(N_PARTICLES) / log(4) + 3)

#define X_MIN -BOUNDS
#define X_MAX BOUNDS
#define Y_MIN -BOUNDS
#define Y_MAX BOUNDS



#include <cuda.h>
#include <cuda_runtime_api.h>

class Particles
{
    public:
        double m_x;
        double m_y;
        double m_vx;
        double m_vy;
        double m_ax;
        double m_ay;
        double m_mass;
};

class QuadTree
{
    public:
        int m_index_node, m_begin_particle, m_end_particle;
        double m_min_x, m_max_x, m_min_y, m_max_y;
};

void initializeParticles(Particles* allParticles, dim3 gridSize, dim3 blockSize);

#endif