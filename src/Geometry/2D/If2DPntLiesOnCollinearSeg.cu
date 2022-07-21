#include "Geometry/2D/If2DPntLiesOnCollinearSeg.cuh"

// ====================================================
// NAME:        If2DPntLiesOnCollinearSeg
// DESCRIPTION: check if a 2D point lies on a collinear segment
//              point q, line segment: p-r
// AUTHOR:      Tingchang YIN
// DATE:        07/05/2022
// ====================================================
template <typename T>
__device__ __host__ bool cuDFNsys::If2DPntLiesOnCollinearSeg(cuDFNsys::Vector2<T> p, cuDFNsys::Vector2<T> q, cuDFNsys::Vector2<T> r)
{
    if (q.x <= max(p.x, r.x) && q.x >= min(p.x, r.x) &&
        q.y <= max(p.y, r.y) && q.y >= min(p.y, r.y))
        return true;

    return false;
}; // If2DPntLiesOnCollinearSeg
template __device__ __host__ bool cuDFNsys::If2DPntLiesOnCollinearSeg<double>(cuDFNsys::Vector2<double> p, cuDFNsys::Vector2<double> q, cuDFNsys::Vector2<double> r);
template __device__ __host__ bool cuDFNsys::If2DPntLiesOnCollinearSeg<float>(cuDFNsys::Vector2<float> p, cuDFNsys::Vector2<float> q, cuDFNsys::Vector2<float> r);