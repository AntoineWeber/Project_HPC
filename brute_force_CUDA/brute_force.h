
#ifndef _BRUTE_FORCE_H_
#define _BRUTE_FORCE_H_

//#define SAVE true

#define N_PARTICLES 2048
#define GRID_SIZE 16
#define BLOCK_SIZE 16

#define GRID_MAX 5
#define GRID_MIN -GRID_MAX
#define R_MAX 7
#define R_OFFSET 3

#define PARTICLE_MASS 100
#define G 6.67408e-11
#define TIMESTEP 100
#define ITERATIONS 100
#define PI 3.14159265

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
std::vector<Position> fillArray(float *pos_x, float *pos_y);
void writeToFile(std::vector<Position> output);

#endif // _BRUTE_FORCE_H_

