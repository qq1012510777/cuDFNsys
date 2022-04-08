///////////////////////////////////////////////////////////////////
// NAME:              Triangle2DOrientation.cuh
//
// PURPOSE:           return orientation of a 2D triangle
//                    true: clockwise
//                    false: counterclockwise
//
// FUNCTIONS/OBJECTS: Triangle2DOrientation
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../../GlobalDef/GlobalDef.cuh"

namespace cuDFNsys
{
__device__ __host__ bool Triangle2DOrientation(float2 a, float2 b, float2 c);
}; // namespace cuDFNsys