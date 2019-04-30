//********************************************************************************
//*
//*  C++/CUDA nbody 
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

#ifndef __BARNESHUTPARTICLESYSTEM_H__
#define __BARNESHUTPARTICLESYSTEM_H__

#include <cuda.h>
#include <cuda_runtime.h>
#include "SimulationParameters.h"


class BarnesHutParticleSystem 
{
	private:
		SimulationParameters parameters;
		int step;
		int numParticles;
		int numNodes;

		float *h_left;
		float *h_right;
		float *h_bottom;
		float *h_top;

		float *h_mass;
		float *h_x;
		float *h_y;
		float *h_vx;
		float *h_vy;
		float *h_ax;
		float *h_ay;

		int *h_child;
		int *h_start;
		int *h_sorted;
		int *h_count;

		float *d_left;
		float *d_right;
		float *d_bottom;
		float *d_top;
		
		float *d_mass;
		float *d_x;
		float *d_y;
		float *d_vx;
		float *d_vy;
		float *d_ax;
		float *d_ay;
		
		int *d_index;
		int *d_child;
		int *d_start;
		int *d_sorted;
		int *d_count;

		int *d_mutex;  //used for locking 

		cudaEvent_t start, stop; // used for timing

		float *h_output;  //host output array for visualization
		float *d_output;  //device output array for visualization

	public:
		BarnesHutParticleSystem(const SimulationParameters p, const int n);
		BarnesHutParticleSystem(const BarnesHutParticleSystem &system);
		BarnesHutParticleSystem& operator=(const BarnesHutParticleSystem &system);
		~BarnesHutParticleSystem();

		int getNumParticles();
		void update();
		void reset();

		const float* getOutputBuffer();

	private:
		void plummerModel(float *mass, float *x, float* y, float *x_vel, float *y_vel, float *x_acc, float *y_acc, int n);
		void diskModel(float *mass, float *x, float* y, float *x_vel, float *y_vel, float *x_acc, float *y_acc, int n);
		void collidingDiskModel(float *mass, float *x, float* y, float *x_vel, float *y_vel, float *x_acc, float *y_acc, int n);
 
};



#endif
