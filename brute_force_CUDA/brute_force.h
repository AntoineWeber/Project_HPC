
#ifndef _BRUTE_FORCE_H_
#define _BRUTE_FORCE_H_

#define N_PARTICLES 32768
#define GRID_SIZE 64
#define BLOCK_SIZE 64

#define GRID_MAX 5
#define GRID_MIN -GRID_MAX

#define PARTICLE_MASS 10
#define G 6.67408e-11
#define TIMESTEP 100
#define ITERATIONS 100

#include <chrono>
#include <iostream>
#include <fstream>

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

void initiateDevice();

#endif // _BRUTE_FORCE_H_

