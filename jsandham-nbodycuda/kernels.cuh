//********************************************************************************
//*
//*  C++/CUDA Barnes Hut
//*  James Sandham
//*  16 May 2016
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

#ifndef __KERNELS_H__
#define __KERNELS_H__

__global__ void set_draw_array_kernel(float *ptr, float *x, float *y, int n);
__global__ void reset_arrays_kernel(int *mutex, float *x, float *y, float *mass, int *count, int *start, int *sorted, int *child, int *index, float *left, float *right, float *bottom, float *top, int n, int m);
__global__ void compute_bounding_box_kernel(int *mutex, float *x, float *y, float *left, float *right, float *bottom, float *top, int n);
__global__ void build_tree_kernel(float *x, float *y, float *mass, int *count, int *start, int *child, int *index, float *left, float *right, float *bottom, float *top, int n, int m);
__global__ void centre_of_mass_kernel(float *x, float *y, float *mass, int *index, int n);
__global__ void sort_kernel(int *count, int *start, int *sorted, int *child, int *index, int n);
__global__ void compute_forces_kernel(float* x, float *y, float *vx, float *vy, float *ax, float *ay, float *mass, int *sorted, int *child, float *left, float *right, int n, float g);
__global__ void update_kernel(float *x, float *y, float *vx, float *vy, float *ax, float *ay, int n, float dt, float d);
__global__ void copy_kernel(float* x, float* y, float* out, int n);

#endif