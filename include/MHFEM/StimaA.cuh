///////////////////////////////////////////////////////////////////
// NAME:              StimaA.cuh
//
// PURPOSE:           Assemble submatrix A for each element
//
// FUNCTIONS/OBJECTS: StimaA
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../Mesh/Mesh.cuh"

namespace cuDFNsys
{
__device__ __host__ void StimaA(cuDFNsys::EleCoor coord, float A[3][3]);
}; // namespace cuDFNsys