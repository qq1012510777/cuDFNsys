///////////////////////////////////////////////////////////////////
// NAME:              OrientationThree2DPnts.cuh
//
// PURPOSE:           return orientation of three 2D points
//                    0: collinear
//                    1: clockwise
//                    2: counterclockwise
//
// FUNCTIONS/OBJECTS: OrientationThree2DPnts
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../../GlobalDef/GlobalDef.cuh"

namespace cuDFNsys
{
__device__ __host__ uint OrientationThree2DPnts(float2 p, float2 q, float2 r, float _tol_);
}; // namespace cuDFNsys