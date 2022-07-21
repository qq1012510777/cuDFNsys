///////////////////////////////////////////////////////////////////
// NAME:              RandomFunction.cuh
//
// PURPOSE:           some functions to generate random numbers
//
// FUNCTIONS/OBJECTS: N/A
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../GlobalDef/GlobalDef.cuh"

namespace cuDFNsys
{
// power law
template <typename T>
__device__ __host__ T RandomPowerlaw(T x0,
                                     T x1,
                                     T alpha_g,
                                     T rand_0_1);

// uniform
template <typename T>
__device__ __host__ T RandomUniform(T a,
                                    T b,
                                    T rand_0_1);

// fisher
template <typename T>
__device__ __host__ T RandomFisher(T rand_0_1,
                                   T fisher_k);

// lognormal: mean is the mean of log(x),
// variance is the variance of log(x)
template <typename T>
__device__ __host__ T RandomLognormal(T mean,
                                      T variance,
                                      T min,
                                      T max,
                                      T rand_0_1);

// return CDF of standard normal distribution (-inf, x)
template <typename T>
__device__ __host__ T StandardNormalCDF(T x);

// return a value that is calculated by inverting the standard normal CDF.
template <typename T>
__device__ __host__ T StandardNormalCDFInv(T p);

// a stable algorithm to calculate polynomial with 8th degrees
template <typename T>
__device__ __host__ T R8PolyValueHorner(int m, T c[], T x);
}; // namespace cuDFNsys