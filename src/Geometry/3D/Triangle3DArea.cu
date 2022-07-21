#include "Geometry/3D/Triangle3DArea.cuh"

// ====================================================
// NAME:        Triangle3DArea
// DESCRIPTION: get 3D triangle area
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================
template <typename T>
__device__ __host__ T cuDFNsys::Triangle3DArea(cuDFNsys::Vector3<T> Pnt1,
                                               cuDFNsys::Vector3<T> Pnt2,
                                               cuDFNsys::Vector3<T> Pnt3)
{
    cuDFNsys::Vector3<T> Pntsets[3] = {Pnt1, Pnt2, Pnt3};

    T area = 0;

    T L[3] = {0};

    for (int i = 0; i < 3; ++i)
    {
        cuDFNsys::Vector3<T> KK = cuDFNsys::MakeVector3<T>(Pntsets[i].x - Pntsets[(i + 1) % 3].x,
                                                           Pntsets[i].y - Pntsets[(i + 1) % 3].y,
                                                           Pntsets[i].z - Pntsets[(i + 1) % 3].z);

        L[i] = pow(KK.x * KK.x + KK.y * KK.y + KK.z * KK.z, 0.5);
    }

    T P = (L[0] + L[1] + L[2]) * 0.5;
    area = pow(P * (P - L[0]) * (P - L[1]) * (P - L[2]), 0.5);

    return area;
}; // Triangle3DArea
template __device__ __host__ double cuDFNsys::Triangle3DArea<double>(cuDFNsys::Vector3<double> Pnt1,
                                                                     cuDFNsys::Vector3<double> Pnt2,
                                                                     cuDFNsys::Vector3<double> Pnt3);
template __device__ __host__ float cuDFNsys::Triangle3DArea<float>(cuDFNsys::Vector3<float> Pnt1,
                                                                   cuDFNsys::Vector3<float> Pnt2,
                                                                   cuDFNsys::Vector3<float> Pnt3);