# Final Project, MATH-454, Parallel and high-performance computing
Author : **Antoine Weber**

The goal of this project was to solve the N-Body problem with the Barnes-Hut algorithm on the GPU using the CUDA platform.

Different milestones were performed before trying the actual resolution of the NBody problem using the Barnes-Hut approximation and the CUDA platform.

* Implementation of the [brute force](brute_force_CUDA) method using the CUDA platform
* Implementation of the [sequential](Sequential_CPU) Barnes-Hut approximation to be ran on the host.
* Implementation of the [Barnes-Hut approximation on the CUDA platform](Barnes_Hut_CUDA).

The parallelization of the Barnes-Hut on the CUDA platform was successfully implemented inducing a linear scaling w.r.t the number of spawned threads. However, the quadtree building was not implemented in a cuda kernel as it was observed that, with the used and designed data structure, the parallelization of the quadtree building would not be worth it in terms of time of work / reward in speedup ratio.

## Compilation
For the [brute force](brute_force_CUDA), a Makefile is furnished for both local compilation and remote using the [EPFL DENEB](https://scitas.epfl.ch/hardware/deneb/) cluster. You may have to fine tune the path to cuda and nvcc within the [Makefile_local](brute_force_CUDA/Makefile_local) if compiling for local use.
```
make -f Makefile_local clean
make -f Makefile_local
```

For the [sequential](Sequential_CPU) code for barnes-hut, a simple CmakeLists.txt is furnished.
```
cmake .
make
```

Finally for the Barnes-Hut in CUDA, the same procedure than with the brute force method can be applied.
