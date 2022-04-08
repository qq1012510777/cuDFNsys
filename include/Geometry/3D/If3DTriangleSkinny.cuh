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

namespace cuDFNsys
{
__device__ __host__ bool If3DTriangleSkinny(float3 A,
                                            float3 B,
                                            float3 C,
                                            float _TOL_s);
}; // namespace cuDFNsys