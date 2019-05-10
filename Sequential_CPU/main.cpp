
#include <chrono>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "barnesHut.hpp"


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


////////////////////////////////////////////////////////////////////////////////
// Program main : 2D
////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
    Timer t1;
    // used for output saving
    #ifdef SAVE
        //open file
        std::ofstream myFile("trajectories.txt");
    #endif

    t1.reset();
    Particles allParticles;
    allParticles.initialize("circle");
    double elapsed_ini = t1.elapsed();

    t1.reset();
    for (unsigned int i=0; i<N_ITERATIONS; i++)
    {
        std::cout << "iteration " << i << std::endl;
        allParticles.resetTree();
        allParticles.buildTree();
        allParticles.computeDisplacement();
        #ifdef SAVE
            allParticles.saveToFile(&myFile);
        #endif
    }
    double elapsed_compute = t1.elapsed();
    allParticles.resetTree();

    std::cout << "The initialisation took : " << elapsed_ini << "seconds." << std::endl;
    std::cout << "While the computation of the " << N_ITERATIONS << " steps took " << elapsed_compute << "seconds." << std::endl;
    std::cout << "with " << N_PARTICULES << " particles." << std::endl;

    return 1;
}


