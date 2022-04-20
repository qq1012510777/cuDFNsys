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
__device__ __host__ uint OrientationThree2DPnts(float2 p, float2 q, float2 r)
{
    int val = (q.y - p.y) * (r.x - q.x) -
              (q.x - p.x) * (r.y - q.y);

    if (val == 0)
        return 0; // collinear

    return (val > 0) ? 1 : 2; // clock or counterclock wise
};
}; // namespace cuDFNsys