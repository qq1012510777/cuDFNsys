///////////////////////////////////////////////////////////////////
// NAME:              MatrixManipulation.cuh
//
// PURPOSE:           Matrix manipulation
//
// FUNCTIONS/OBJECTS: N/A
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../GlobalDef/GlobalDef.cuh"

namespace cuDFNsys
{
// cross product of two float3 variables
__device__ __host__ float3 CrossProductFloat3(float3 v1, float3 v2);
// cross product of two float2 variables
__device__ __host__ float CrossProductFloat2(float2 A, float2 B);
// product of a square matrix and a column vector
__device__ __host__ float3 ProductSquare3Float3(float A[3][3], float3 B);
}; // namespace cuDFNsys