
#include "call_kernels.cuh"
#include "kernels.cuh"

__global__ void initialize_particles_circle(Particles* allParticles)
{
    int i = threadIdx.x + blockIdx.x*blockDim.x;
    int stride = gridDim.x * blockDim.x;
    int offset = 0;
    
    // to initialize the cuda rand
    curandState state;
    curand_init(clock64(), i, 0, &state);

    while (i + offset < N_PARTICLES)
    {
        allParticles[i + offset].m_mass = (double)MASS;

        double r = CIRCLE_OFFSET + curand_uniform(&state)*2*BOUNDS;
        double alpha = curand_uniform(&state)*2*M_PI;

        allParticles[i + offset].m_x = r*cos(alpha);
        allParticles[i + offset].m_y = r*sin(alpha);

        // set velocity to 0
        allParticles[i + offset].m_vx = 0.0;
        allParticles[i + offset].m_vy = 0.0;

        // set acceleration to 0
        allParticles[i + offset].m_ax = 0.0;
        allParticles[i + offset].m_ay = 0.0;

        offset += stride;
    }
}

__global__ void build_tree(Particles* allParticles, QuadTree* allNodes, int depth, int n_level_nodes, int blockind)
{
    printf("block : %d \n", blockIdx.x);
    int identifier = depth*(blockIdx.x)*n_level_nodes + blockIdx.x;

    if (depth == 3)
    {
        // fill what needs to be filled
        //printf("max depth \n");
        return;
    }

    //printf("identifier : %d \n", identifier);
    build_tree<<<4,1>>>(allParticles, allNodes, depth+1, n_level_nodes*4);
}