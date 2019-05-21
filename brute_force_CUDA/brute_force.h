
#ifndef _BRUTE_FORCE_H_
#define _BRUTE_FORCE_H_

//#define SAVE true

#define N_PARTICLES 16384
#define GRID_SIZE 256
#define BLOCK_SIZE 256

#define GRID_MAX 5
#define GRID_MIN -GRID_MAX
#define R_MAX 7
#define R_OFFSET 3

#define G 1.0
#define PARTICLE_MASS 0.001
#define TIMESTEP 0.001
#define ITERATIONS 100
#define EPSILON 0.001

#define SIMU_BOUND_X 20
#define SIMU_BOUND_Y 20

#include <chrono>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

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

struct Position
{
    double x,y;
};

void initiateDevice();

#ifdef SAVE
std::vector<Position> fillArray(double *pos_x, double *pos_y);
#endif

#endif // _BRUTE_FORCE_H_

