#include "Geometry/2D/DistancePnt2DSeg.cuh"

// ====================================================
// NAME:        DistancePnt2DSeg
// DESCRIPTION: get distance between 2D point and 
//              2D segment
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================
__device__ __host__ float cuDFNsys::DistancePnt2DSeg(float2 pnt,
                                                     float2 *verts)
{
    float2 AB = make_float2(verts[1].x - verts[0].x, verts[1].y - verts[0].y);

    float2 AE = make_float2(pnt.x - verts[0].x, pnt.y - verts[0].y);

    float2 BE = make_float2(pnt.x - verts[1].x, pnt.y - verts[1].y);

    float AB_BE = AB.x * BE.x + AB.y * BE.y;
    float AB_AE = AB.x * AE.x + AB.y * AE.y;

    float reqAns = 0;

    if (AB_BE > 0)
        reqAns = sqrt(BE.x * BE.x + BE.y * BE.y);
    else if (AB_AE < 0)
        reqAns = sqrt(AE.x * AE.x + AE.y * AE.y);
    else
    {
        float x1 = AB.x;
        float y1 = AB.y;
        float x2 = AE.x;
        float y2 = AE.y;
        float mod = sqrt(x1 * x1 + y1 * y1);
        reqAns = abs(x1 * y2 - y1 * x2) / mod;
    };

    return reqAns;
};