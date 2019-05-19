

#include "call_kernels.cuh"
#include "kernels.cuh"



void initializeParticles(Particles* allParticles, dim3 gridSize, dim3 blockSize)
{
	initialize_particles_circle<<<gridSize, blockSize>>>(allParticles);
}

void constructTree(Particles* allParticles,QuadTree* allNodes, int depth)
{
	build_tree<<<1,128>>>(allParticles, allNodes, depth);
}