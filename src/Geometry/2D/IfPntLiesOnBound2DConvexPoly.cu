#include "Geometry/2D/IfPntLiesOnBound2DConvexPoly.cuh"

// ====================================================
// NAME:        IfPntLiesOnBound2DConvexPoly
// DESCRIPTION: if a point lies on 2D bound
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================
__device__ __host__ bool cuDFNsys::IfPntLiesOnBound2DConvexPoly(float2 pnt,
                                                                float2 *verts,
                                                                int N,
                                                                float _tol_)
{
    for (int i = 0; i < N; ++i)
    {
        float2 Seg[2] = {verts[i], verts[(i + 1) % N]};

        float dist = cuDFNsys::DistancePnt2DSeg(pnt, Seg);

        if (dist < _tol_)
        {
            return true;
        }
    }

    return false;
}; // IfPntLiesOnBound2DConvexPoly

// ====================================================
// NAME:        IfPntLiesOnBound2DConvexPolyReturnEdgeNO
// DESCRIPTION: if a point lies on 2D bound
//              also return the edge NO
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================
__device__ __host__ bool cuDFNsys::IfPntLiesOnBound2DConvexPolyReturnEdgeNO(float2 pnt,
                                                                            float2 *verts,
                                                                            int N,
                                                                            float _tol_,
                                                                            int *edgeNO) // 0 1 2
{
    for (int i = 0; i < N; ++i)
    {
        float2 Seg[2] = {verts[i], verts[(i + 1) % N]};

        float dist = cuDFNsys::DistancePnt2DSeg(pnt, Seg);

        if (dist < _tol_)
        {
            *edgeNO = i;
            return true;
        }
    }

    return false;
}; // IfPntLiesOnBound2DConvexPolyReturnEdgeNO