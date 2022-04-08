#include "Geometry/3D/If3DTriangleSkinny.cuh"

// ====================================================
// NAME:        If3DTriangleSkinny
// DESCRIPTION: If a triangle is skinny
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================
__device__ __host__ bool cuDFNsys::If3DTriangleSkinny(float3 A,
                                                      float3 B,
                                                      float3 C,
                                                      float _TOL_s)
{
    float3 Tri[3] = {A, B, C};

    for (int i = 0; i < 3; ++i)
    {
        float x1, x2, y1, y2, z1, z2;
        x1 = Tri[(i + 1) % 3].x - Tri[i].x;
        y1 = Tri[(i + 1) % 3].y - Tri[i].y;
        z1 = Tri[(i + 1) % 3].z - Tri[i].z;

        x2 = Tri[(i + 2) % 3].x - Tri[(i + 1) % 3].x;
        y2 = Tri[(i + 2) % 3].y - Tri[(i + 1) % 3].y;
        z2 = Tri[(i + 2) % 3].z - Tri[(i + 1) % 3].z;

        float dot = x1 * x2 + y1 * y2 + z1 * z2; //    #between [x1, y1, z1] and [x2, y2, z2]
        float lenSq1 = x1 * x1 + y1 * y1 + z1 * z1;
        float lenSq2 = x2 * x2 + y2 * y2 + z2 * z2;      //
        float angle = acos(dot / sqrt(lenSq1 * lenSq2)); //

        angle = angle * 180.0f / M_PI;

        if (abs(angle - 0) < _TOL_s ||
            abs(angle - 180.0f) < _TOL_s)
            return true;
    }
    return false;
}; // If3DTriangleSkinny
