///////////////////////////////////////////////////////////////////
// NAME:              If3DTriangleSkinny.cuh
//
// PURPOSE:           if a triangle is skinny?
//
// FUNCTIONS/OBJECTS: If3DTriangleSkinny
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../../GlobalDef/GlobalDef.cuh"
#include "./DataTypeSelector/DataTypeSelector.cuh"

namespace cuDFNsys
{
template <typename T>
__device__ __host__ bool If3DTriangleSkinny(cuDFNsys::Vector3<T> A,
                                            cuDFNsys::Vector3<T> B,
                                            cuDFNsys::Vector3<T> C,
                                            T _TOL_s);
}; // namespace cuDFNsys