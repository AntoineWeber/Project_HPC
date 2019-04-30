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
#ifndef __NBODYVISUALIZER_H__
#define __NBODYVISUALIZER_H__


#include <iostream>
#include "SimulationParameters.h"
#include "BarnesHutParticleSystem.h"


#define GLEW_STATIC
#include <GL/glew.h>
#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>

#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "glm/gtc/type_ptr.hpp"



class NBodyVisualizer
{
	private:
		int numOfBodies;
		BarnesHutParticleSystem *particles;
		SimulationParameters parameters;

		sf::ContextSettings *settings;
		sf::Window *window;

		GLuint vao;
		GLuint vbo;

		GLuint vertexShader;
		GLuint fragmentShader;
		GLuint shaderProgram;

		void displayDeviceProperties();

	public:
		NBodyVisualizer(const SimulationParameters p, const int numBodies);
		NBodyVisualizer(const NBodyVisualizer &visualizer);
		NBodyVisualizer& operator=(const NBodyVisualizer &visualizer);
		~NBodyVisualizer();

		void runSimulation();
};



#endif