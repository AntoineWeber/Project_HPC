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


#ifndef __DEBUG_H__
#define __DEBUG_H__


#define DEBUG 1
#define TIMING 0

#define DEBUG_PRINT(stream, statement) \
	do { if(DEBUG) (stream) << "DEBUG: "<< __FILE__<<"("<<__LINE__<<") " << (statement) << std::endl;} while(0)

#define TIMING_PRINT(stream, statement, time) \
	do { if(TIMING) (stream) << "TIMING: "<<__FILE__<<"("<<__LINE__<<") " << (statement) << (time) << std::endl;} while(0)


void DEBUG_RUN_TESTS(float *x, float *y, float *mass, int *count, int *start, int *sorted, int *child, float *left, float *right, float *bottom, float *top, int n, int m);

#endif