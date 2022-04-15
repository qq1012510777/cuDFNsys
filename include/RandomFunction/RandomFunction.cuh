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
__device__ __host__ float RandomPowerlaw(float x0,
                                         float x1,
                                         float alpha_g,
                                         float rand_0_1);
// uniform
__device__ __host__ float RandomUniform(float a,
                                        float b,
                                        float rand_0_1);
// fisher
__device__ __host__ float RandomFisher(float rand_0_1,
                                       float fisher_k);
// lognormal: mean is the mean of log(x),
// variance is the variance of log(x)
__device__ __host__ float RandomLognormal(float mean,
                                          float variance,
                                          float min,
                                          float max,
                                          float rand_0_1);
// return CDF of standard normal distribution (-inf, x)
__device__ __host__ float StandardNormalCDF(float x);
// return a value that is calculated by inverting the standard normal CDF.
__device__ __host__ float StandardNormalCDFInv(float p);
// a stable algorithm to calculate polynomial with 8th degrees
__device__ __host__ float R8PolyValueHorner(int m, float c[], float x);
}; // namespace cuDFNsys