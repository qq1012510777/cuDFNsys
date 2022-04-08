///////////////////////////////////////////////////////////////////
// NAME:              DistancePnt2DSeg.cuh
//
// PURPOSE:           get distance between 2D point and 2D segment
//
// FUNCTIONS/OBJECTS: DistancePnt2DSeg
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../../GlobalDef/GlobalDef.cuh"

namespace cuDFNsys
{
__device__ __host__ float DistancePnt2DSeg(float2 pnt,
                                           float2 *verts);
}; // namespace cuDFNsys