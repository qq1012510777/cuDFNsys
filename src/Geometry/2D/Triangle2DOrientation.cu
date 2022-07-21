#include "Geometry/2D/Triangle2DOrientation.cuh"

// ====================================================
// NAME:        Triangle2DOrientation
// DESCRIPTION: return orientation of a 2D triangle
//              true: clockwise
//              false: counterclockwise
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================
template <typename T>
__device__ __host__ bool cuDFNsys::Triangle2DOrientation(cuDFNsys::Vector2<T> a,
                                                         cuDFNsys::Vector2<T> b,
                                                         cuDFNsys::Vector2<T> c)
{
    T d = (b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y);

    if (d < 0)
        return true;
    else
        return false;
}; // Triangle2DOrientation
template __device__ __host__ bool cuDFNsys::Triangle2DOrientation<double>(cuDFNsys::Vector2<double> a,
                                                                          cuDFNsys::Vector2<double> b,
                                                                          cuDFNsys::Vector2<double> c);
template __device__ __host__ bool cuDFNsys::Triangle2DOrientation<float>(cuDFNsys::Vector2<float> a,
                                                                         cuDFNsys::Vector2<float> b,
                                                                         cuDFNsys::Vector2<float> c);