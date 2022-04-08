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
}; // namespace cuDFNsys