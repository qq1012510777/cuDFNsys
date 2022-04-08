///////////////////////////////////////////////////////////////////
// NAME:              Intersection2DLine2DSeg.cuh
//
// PURPOSE:           Identify intersection between
//                    2D line (infinite) and 2D segment
//
// FUNCTIONS/OBJECTS: Intersection2DLine2DSeg
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////

#pragma once
#include "../../GlobalDef/GlobalDef.cuh"

namespace cuDFNsys
{
// Identify intersection between 2D line (infinite) and 2D segment
__device__ __host__ bool Intersection2DLine2DSeg(float2 *Line,
                                                 float2 *Seg,
                                                 int *sign_, // 1, pnt; 2, seg; 3, none
                                                 float2 *intersection,
                                                 float _TOL_);
}; // namespace cuDFNsys
