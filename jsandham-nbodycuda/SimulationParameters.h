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

#ifndef __SIMULATIONPARAMETERS_H__
#define __SIMULATIONPARAMETERS_H__



typedef enum Model
{
	disk_model,
	colliding_disk_model,
	plummer_model
}Model;


typedef struct SimulationParameters
{
	Model model;
	bool opengl;
	bool debug;
	bool benchmark;
	bool fullscreen;
	float iterations;
	float timestep;
	float gravity;
	float dampening;
}SimulationParameters;

#endif
