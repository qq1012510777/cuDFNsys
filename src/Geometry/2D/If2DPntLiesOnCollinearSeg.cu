#include "Geometry/2D/If2DPntLiesOnCollinearSeg.cuh"

// ====================================================
// NAME:        If2DPntLiesOnCollinearSeg
// DESCRIPTION: check if a 2D point lies on a collinear segment
//              point q, line segment: p-r
// AUTHOR:      Tingchang YIN
// DATE:        07/05/2022
// ====================================================

__device__ __host__ bool cuDFNsys::If2DPntLiesOnCollinearSeg(float2 p, float2 q, float2 r)
{
    if (q.x <= max(p.x, r.x) && q.x >= min(p.x, r.x) &&
        q.y <= max(p.y, r.y) && q.y >= min(p.y, r.y))
        return true;

    return false;
}; // If2DPntLiesOnCollinearSeg