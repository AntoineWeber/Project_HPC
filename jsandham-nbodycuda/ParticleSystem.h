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
// NBodyCuda is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.
//
//********************************************************************************

#ifndef __PARTICLESYSTEM_H__
#define __PARTICLESYSTEM_H__

#include <iostream>
#include "SimulationParameters.h"

class ParticleSystem
{
	public:
		SimulationParameters parameters;
		
		ParticleSystem(const SimulationParameters p, const int n){}
		virtual ~ParticleSystem(){};

		virtual int getNumParticles() = 0; 
		virtual void update() = 0;
		virtual void reset() = 0;
		virtual float* getOutputBuffer() = 0;

		void plummerModel(float *mass, float *x, float* y, float *x_vel, float *y_vel, float *x_acc, float *y_acc, int n);
		void diskModel(float *mass, float *x, float* y, float *x_vel, float *y_vel, float *x_acc, float *y_acc, int n);
		void collidingDiskModel(float *mass, float *x, float* y, float *x_vel, float *y_vel, float *x_acc, float *y_acc, int n);
};


#endif