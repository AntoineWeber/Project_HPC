
#include <iostream>
#include <fstream>
#include <vector>

#include "barnesHut.hpp"


////////////////////////////////////////////////////////////////////////////////
// Program main : 2D
////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
    Particles allParticles;
    allParticles.initialize("uniform");

    for (unsigned int i=0; i<N_ITERATIONS; i++)
    {
        std::cout << "iteration " << i << std::endl;
        allParticles.resetTree();
        allParticles.computeBoundingBox();
        allParticles.buildTree();
        allParticles.computeForce();
    }

    allParticles.resetTree();

    return 1;
}


