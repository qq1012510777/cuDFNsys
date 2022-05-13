#include "Geometry/2D/OrientationThree2DPnts.cuh"
// ====================================================
// NAME:        OrientationThree2DPnts
// DESCRIPTION: return orientation of three 2D points
//                    0: collinear
//                    1: clockwise
//                    2: counterclockwise
// AUTHOR:      Tingchang YIN
// DATE:        07/05/2022
// ====================================================

__device__ __host__ uint cuDFNsys::OrientationThree2DPnts(float2 p, float2 q, float2 r, float _tol_)
{
    float val = (q.y - p.y) * (r.x - q.x) -
                (q.x - p.x) * (r.y - q.y);

    if (abs(val) < _tol_)
        return 0; // collinear

    // clock or counterclock wise
    return (val > 0) ? 1 : 2;
}; // OrientationThree2DPnts