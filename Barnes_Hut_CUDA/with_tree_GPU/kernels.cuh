
#ifndef __KERNELS_H__
#define __KERNELS_H__

#include <stdio.h>
#include <math.h>
#include <curand.h>
#include <curand_kernel.h>
#include <cuda.h>


#include "call_kernels.cuh"

__global__ void initialize_particles_circle(Particles* allParticles);

#endif