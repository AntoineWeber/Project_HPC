

#include "call_kernels.cuh"
#include "kernels.cuh"



void initializeParticles(Particles* allParticles, dim3 gridSize, dim3 blockSize)
{
	initialize_particles_circle<<<gridSize, blockSize>>>(allParticles);
}

void constructTree(Particles* allParticles,QuadTree* allNodes, int depth, int n_level_nodes, int blockind)
{
	build_tree<<<1,1>>>(allParticles, allNodes, depth, n_level_nodes, blockind);
}