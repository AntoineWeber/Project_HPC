#Makefile for nbody code

#define variables
objects= main.o NBodyVisualizer.o kernels.o BarnesHutParticleSystem.o Particle_cuda.o debug.o
NVCC= nvcc              #cuda c compiler
CPP= g++                 #c++ compiler
opt= -O2 -g -G           #optimization flag
ARCH= -arch=sm_61        #cuda compute capability
#LIBS= -lglut -lGLU -lGL  
#LIBS= -lGLU -lGL -lGLEW -lm -lsfml-graphics -lsfml-window -lsfml-system 
LIBS= -lGL -lGLEW -lsfml-graphics -lsfml-window -lsfml-system
INCLUDE = -I/home/antoine/glm-0.9.9.5/glm/
execname= main

build: $(objects)
	$(NVCC) $(opt) -o $(execname) $(objects) $(LIBS) 


kernels.o: kernels.cu
	$(NVCC) $(opt) $(ARCH) -maxrregcount=32 -c kernels.cu
BarnesHutParticleSystem.o: BarnesHutParticleSystem.cpp
	$(NVCC) $(opt) -std=c++11 -c BarnesHutParticleSystem.cpp 
Particle_cuda.o: Particle_cuda.cu
	$(NVCC) $(opt) $(ARCH) -c Particle_cuda.cu 
debug.o: debug.cpp
	$(NVCC) $(opt) $(ARCH) -c debug.cpp
NBodyVisualizer.o: NBodyVisualizer.cpp
	$(NVCC) $(opt) $(ARCH) -c $(INCLUDE) NBodyVisualizer.cpp 
main.o: main.cpp 
	$(NVCC) $(opt) $(ARCH) -c $(INCLUDE) main.cpp


clean:
	rm $(objects)
