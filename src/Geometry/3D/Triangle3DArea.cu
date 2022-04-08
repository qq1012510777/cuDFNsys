#include "Geometry/3D/Triangle3DArea.cuh"

// ====================================================
// NAME:        Triangle3DArea
// DESCRIPTION: get 3D triangle area
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================

__device__ __host__ float cuDFNsys::Triangle3DArea(float3 Pnt1,
                                                   float3 Pnt2,
                                                   float3 Pnt3)
{
    float3 Pntsets[3] = {Pnt1, Pnt2, Pnt3};

    float area = 0;

    float L[3] = {0};

    for (int i = 0; i < 3; ++i)
    {
        float3 KK = make_float3(Pntsets[i].x - Pntsets[(i + 1) % 3].x,
                                Pntsets[i].y - Pntsets[(i + 1) % 3].y,
                                Pntsets[i].z - Pntsets[(i + 1) % 3].z);

        L[i] = pow(KK.x * KK.x + KK.y * KK.y + KK.z * KK.z, 0.5);
    }

    float P = (L[0] + L[1] + L[2]) * 0.5;
    area = pow(P * (P - L[0]) * (P - L[1]) * (P - L[2]), 0.5);

    return area;
}; // Triangle3DArea