///////////////////////////////////////////////////////////////////
// NAME:              Triangle2DArea.cuh
//
// PURPOSE:           Get 2D triangle area
//
// FUNCTIONS/OBJECTS: Triangle2DArea
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../../GlobalDef/GlobalDef.cuh"

namespace cuDFNsys
{
__device__ __host__ float Triangle2DArea(float2 Pnt1,
                                         float2 Pnt2,
                                         float2 Pnt3);
}; // namespace cuDFNsys