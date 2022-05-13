///////////////////////////////////////////////////////////////////
// NAME:              If2DPntLiesOnCollinearSeg.cuh
//
// PURPOSE:           check if a 2D point lies on a collinear segment
//                    point q, line segment: p-r
//
// FUNCTIONS/OBJECTS: If2DPntLiesOnCollinearSeg
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../../GlobalDef/GlobalDef.cuh"

namespace cuDFNsys
{
__device__ __host__ bool If2DPntLiesOnCollinearSeg(float2 p, float2 q, float2 r);
}; // namespace cuDFNsys