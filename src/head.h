#pragma once

#include <omp.h>
#include <stdlib.h>

// #include <stdexcept>
// #define Thread_Num 16
// // if define Thread_Nums in makefile
// // e.g., CXXFLAGS = -DThread_Nums=4 -fopenmp

// you can define Thread_Num in makefile:
// FLAGS += -DThread_Num=16 -fopenmp
#ifndef Thread_Num
// tells the compiler Thread_Num variable is defined somewhere else
// if used in code, need to define it
// you can define Thread_Num in .cpp file
// then use .ini file to store Thread_Num
extern int Thread_Num;
#endif


// float number type
// FLAGS        += -DFLOAT_TYPE=double
#ifdef FLOAT_TYPE
using FLOAT = FLOAT_TYPE;
#else
using FLOAT = double;
#endif

// platform
#ifdef __CUDACC__

#define __any__ __device__ __host__

#ifdef __CUDA_ARCH__

#define __constexpr__ __constant__ constexpr 

#else

#define __constexpr__ constexpr

#endif

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/tuple.h>


#else

#define __any__
#define __constexpr__ constexpr
#define __global__
#define __host__
#define __device__


#endif


// users can decide whether to add the debug macros:
// FLAGS        += -D_FLB_DEBUG