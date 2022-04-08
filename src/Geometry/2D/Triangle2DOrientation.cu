#include "Geometry/2D/Triangle2DOrientation.cuh"

// ====================================================
// NAME:        Triangle2DOrientation
// DESCRIPTION: return orientation of a 2D triangle
//              true: clockwise
//              false: counterclockwise
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================
__device__ __host__ bool cuDFNsys::Triangle2DOrientation(float2 a, float2 b, float2 c)
{
    float d = (b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y);

    if (d < 0)
        return true;
    else
        return false;
}; // Triangle2DOrientation