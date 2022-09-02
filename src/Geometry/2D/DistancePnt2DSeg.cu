#include "Geometry/2D/DistancePnt2DSeg.cuh"

// ====================================================
// NAME:        DistancePnt2DSeg
// DESCRIPTION: get distance between 2D point and
//              2D segment
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================
template <typename T>
__device__ __host__ T cuDFNsys::DistancePnt2DSeg(cuDFNsys::Vector2<T> pnt, cuDFNsys::Vector2<T> *verts)
{
    // printf("pnt: %.40f, %.40f\n", pnt.x, pnt.y);
    // printf("verts[0]: %.40f, %.40f\n", verts[0].x, verts[0].y);
    // printf("verts[1]: %.40f, %.40f\n", verts[1].x, verts[1].y);
    cuDFNsys::Vector2<T> AB = cuDFNsys::MakeVector2(verts[1].x - verts[0].x, verts[1].y - verts[0].y);

    cuDFNsys::Vector2<T> AE = cuDFNsys::MakeVector2(pnt.x - verts[0].x, pnt.y - verts[0].y);

    cuDFNsys::Vector2<T> BE = cuDFNsys::MakeVector2(pnt.x - verts[1].x, pnt.y - verts[1].y);

    T AB_BE = AB.x * BE.x + AB.y * BE.y;
    T AB_AE = AB.x * AE.x + AB.y * AE.y;
    /// printf("AB_BE: %.40f, AB_AE: %.40f\n", AB_BE, AB_AE);
    T reqAns = 0;

    if (AB_BE > 0)
        reqAns = sqrt(BE.x * BE.x + BE.y * BE.y);
    else if (AB_AE < 0)
        reqAns = sqrt(AE.x * AE.x + AE.y * AE.y);
    else
    {
        T x1 = AB.x;
        T y1 = AB.y;
        T x2 = AE.x;
        T y2 = AE.y;
        T mod = sqrt(x1 * x1 + y1 * y1);
        reqAns = abs(x1 * y2 - y1 * x2) / mod;
    };

    return reqAns;
}; // DistancePnt2DSeg
template __device__ __host__ double cuDFNsys::DistancePnt2DSeg<double>(cuDFNsys::Vector2<double> pnt, cuDFNsys::Vector2<double> *verts);
template __device__ __host__ float cuDFNsys::DistancePnt2DSeg<float>(cuDFNsys::Vector2<float> pnt, cuDFNsys::Vector2<float> *verts);