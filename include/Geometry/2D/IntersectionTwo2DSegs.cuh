///////////////////////////////////////////////////////////////////
// NAME:              IntersectionTwo2DSegs.cuh
//
// PURPOSE:           Identify intersection between two 2D segments
//
// FUNCTIONS/OBJECTS: IntersectionTwo2DSegs
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../../GlobalDef/GlobalDef.cuh"
#include "If2DPntLiesOnCollinearSeg.cuh"
#include "IntersectionTwoCollinearSegs.cuh"
#include "OrientationThree2DPnts.cuh"

namespace cuDFNsys
{
__device__ __host__ bool IntersectionTwo2DSegs(float2 *Seg_1,
                                               float2 *Seg_2,
                                               float2 *intersection,
                                               int *sign, // 1, pnt; 2, seg; 3, none
                                               float _TOL_);
}; // namespace cuDFNsys