
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
        allParticles.resetTree();
        allParticles.computeBoundingBox();
        allParticles.buildTree();
    }

    return 1;
}


