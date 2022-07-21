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
#include "../DataTypeSelector/DataTypeSelector.cuh"

namespace cuDFNsys
{
template <typename T>
__device__ __host__ void StimaA(cuDFNsys::EleCoor<T> coord, T A[3][3]);
}; // namespace cuDFNsys