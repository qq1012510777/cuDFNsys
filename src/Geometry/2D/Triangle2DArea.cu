#include "Geometry/2D/Triangle2DArea.cuh"

// ====================================================
// NAME:        Triangle2DArea
// DESCRIPTION: get 2D triangle area
// AUTHOR:      Tingchang YIN
// DATE:        19/04/2022
// ====================================================

__device__ __host__ float cuDFNsys::Triangle2DArea(float2 Pnt1,
                                                   float2 Pnt2,
                                                   float2 Pnt3)
{
    float2 Pntsets[3] = {Pnt1, Pnt2, Pnt3};

    float area = 0;

    float L[3] = {0};

    for (int i = 0; i < 3; ++i)
    {
        float2 KK = make_float2(Pntsets[i].x - Pntsets[(i + 1) % 3].x,
                                Pntsets[i].y - Pntsets[(i + 1) % 3].y);

        L[i] = pow(KK.x * KK.x + KK.y * KK.y, 0.5);
    }

    float P = (L[0] + L[1] + L[2]) * 0.5;
    area = pow(P * (P - L[0]) * (P - L[1]) * (P - L[2]), 0.5);

    return area;
}; // Triangle2DArea