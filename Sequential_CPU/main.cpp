
#include <chrono>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "barnesHut.hpp"


////////////////////////////////////////////////////////////////////////////////
// Program main : 2D
////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
    // used for output saving
    #ifdef SAVE
        //open file
        std::ofstream myFile("trajectories.txt");
    #endif

    Particles allParticles;
    allParticles.initialize("uniform");

    for (unsigned int i=0; i<N_ITERATIONS; i++)
    {
        std::cout << "iteration " << i << std::endl;
        allParticles.resetTree();
        allParticles.computeBoundingBox();
        allParticles.buildTree();
        allParticles.computeDisplacement();
        #ifdef SAVE
            allParticles.saveToFile(&myFile);
        #endif
    }
    allParticles.resetTree();

    return 1;
}


