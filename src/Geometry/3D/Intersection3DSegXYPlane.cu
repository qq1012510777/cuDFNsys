#include "Geometry/3D/Intersection3DSegXYPlane.cuh"

// ====================================================
// NAME:        Intersection3DSegXYPlane
// DESCRIPTION: Identify intersection between
//              a 3D segment and the XY plane
// AUTHOR:      Tingchang YIN
// DATE:        07/04/2022
// ====================================================

__device__ __host__ bool cuDFNsys::Intersection3DSegXYPlane(float3 *Seg,
                                                            float3 *Intersec_PNT,
                                                            int *sign_,
                                                            float _TOL_) // sign: 1: pnt; 2: seg; -1: none;
{
    if (abs(Seg[0].z) < _TOL_ && abs(Seg[1].z) < _TOL_)
    {
        *sign_ = 2;
        Intersec_PNT[0] = Seg[0];
        Intersec_PNT[1] = Seg[1];
        return true;
    }

    if (abs(Seg[0].z) < _TOL_)
    {
        Intersec_PNT[0] = Seg[0];
        *sign_ = 1;
        return true;
    }
    else if (abs(Seg[1].z) < _TOL_)
    {
        Intersec_PNT[0] = Seg[1];
        *sign_ = 1;
        return true;
    }

    float max_z = 0;
    float min_z = 0;

    if (Seg[0].z <= Seg[1].z)
    {
        max_z = Seg[1].z;
        min_z = Seg[0].z;
    }
    else
    {
        min_z = Seg[1].z;
        max_z = Seg[0].z;
    }

    // common cases: cross the xy plane

    if (max_z > 0 && min_z < 0)
    {
        // cross the xy plane
        float3 center_seg = make_float3(0.5 * (Seg[0].x + Seg[1].x),
                                        0.5 * (Seg[0].y + Seg[1].y),
                                        0.5 * (Seg[0].z + Seg[1].z));

        float3 vector_seg = make_float3((Seg[0].x - Seg[1].x),
                                        (Seg[0].y - Seg[1].y),
                                        (Seg[0].z - Seg[1].z));

        float vpt = vector_seg.z * 1;

        float t = (0 - center_seg.z) * 1 / vpt;

        float3 Pnt;
        Pnt.x = center_seg.x + vector_seg.x * t;
        Pnt.y = center_seg.y + vector_seg.y * t;
        Pnt.z = center_seg.z + vector_seg.z * t;

        Intersec_PNT[0] = Pnt;

        *sign_ = 1;
        return true;
    }

    *sign_ = -1;
    return false;
}; // Intersection3DSegXYPlane
