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

template <typename T>
__device__ __host__ uint cuDFNsys::OrientationThree2DPnts(cuDFNsys::Vector2<T> p, cuDFNsys::Vector2<T> q, cuDFNsys::Vector2<T> r, T _tol_)
{
    T val = (q.y - p.y) * (r.x - q.x) -
            (q.x - p.x) * (r.y - q.y);

    if (abs(val) < _tol_)
        return 0; // collinear

    // clock or counterclock wise
    return (val > 0) ? 1 : 2;
}; // OrientationThree2DPnts
template __device__ __host__ uint cuDFNsys::OrientationThree2DPnts<double>(cuDFNsys::Vector2<double> p, cuDFNsys::Vector2<double> q, cuDFNsys::Vector2<double> r, double _tol_);
template __device__ __host__ uint cuDFNsys::OrientationThree2DPnts<float>(cuDFNsys::Vector2<float> p, cuDFNsys::Vector2<float> q, cuDFNsys::Vector2<float> r, float _tol_);