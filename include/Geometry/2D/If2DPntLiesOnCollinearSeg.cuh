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
__device__ __host__ bool If2DPntLiesOnCollinearSeg(float2 p, float2 q, float2 r)
{
    if (q.x <= max(p.x, r.x) && q.x >= min(p.x, r.x) &&
        q.y <= max(p.y, r.y) && q.y >= min(p.y, r.y))
        return true;

    return false;
};
}; // namespace cuDFNsys