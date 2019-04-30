//********************************************************************************
//*
//*  C++/CUDA Barnes Hut
//*  James Sandham
//*  128 May 20128
//*
//********************************************************************************

//********************************************************************************
// This software is based loosely on the opencl implementation by Benjamin Neukom 
// and the paper "An Efficient CUDA Implementation of the Tree-Based Barnes Hut 
// n-Body Algorithm" by Martin Burtscher and Keshav Pingali.
//
// Barnes-Hut is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.
//
//********************************************************************************


#include "Particle_cuda.cuh"
#include "kernels.cuh"


dim3 blockSize = 256;
dim3 gridSize = 256;

void SetDrawArray(float *ptr, float *x, float *y, int n)
{
	set_draw_array_kernel<<< gridSize, blockSize>>>(ptr, x, y, n);
}


void ResetArrays(int *mutex, float *x, float *y, float *mass, int *count, int *start, int *sorted, int *child, int *index, float *left, float *right, float *bottom, float *top, int n, int m)
{
	reset_arrays_kernel<<< gridSize, blockSize >>>(mutex, x, y, mass, count, start, sorted, child, index, left, right, bottom, top, n, m);
}


void ComputeBoundingBox(int *mutex, float *x, float *y, float *left, float *right, float *bottom, float *top, int n)
{
	compute_bounding_box_kernel<<< gridSize, blockSize >>>(mutex, x, y, left, right, bottom, top, n);
}


void BuildQuadTree(float *x, float *y, float *mass, int *count, int *start, int *child, int *index, float *left, float *right, float *bottom, float *top, int n, int m)
{
	build_tree_kernel<<< gridSize, blockSize >>>(x, y, mass, count, start, child, index, left, right, bottom, top, n, m);
}


void ComputeCentreOfMass(float *x, float *y, float *mass, int *index, int n)
{
	centre_of_mass_kernel<<<gridSize, blockSize>>>(x, y, mass, index, n);
}


void SortParticles(int *count, int *start, int *sorted, int *child, int *index, int n)
{
	sort_kernel<<< gridSize, blockSize >>>(count, start, sorted, child, index, n);
}


void CalculateForces(float* x, float *y, float *vx, float *vy, float *ax, float *ay, float *mass, int *sorted, int *child, float *left, float *right, int n, float g)
{
	compute_forces_kernel<<< gridSize, blockSize >>>(x, y, vx, vy, ax, ay, mass, sorted, child, left, right, n, g);
}


void IntegrateParticles(float *x, float *y, float *vx, float *vy, float *ax, float *ay, int n, float dt, float d)
{
	update_kernel<<<gridSize, blockSize >>>(x, y, vx, vy, ax, ay, n, dt, d);
}


void FillOutputArray(float *x, float *y, float *out, int n)
{
	copy_kernel<<<gridSize, blockSize >>>(x, y, out, n);
}