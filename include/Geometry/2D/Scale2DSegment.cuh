///////////////////////////////////////////////////////////////////
// NAME:              Scale2DSegment.cuh
//
// PURPOSE:           scale a 2D line segment along the center
//
// FUNCTIONS/OBJECTS: Scale2DSegment
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../../GlobalDef/GlobalDef.cuh"

namespace cuDFNsys
{
__host__ __device__ void Scale2DSegment(float2 *Segment, float scaleF);
}; // namespace cuDFNsys