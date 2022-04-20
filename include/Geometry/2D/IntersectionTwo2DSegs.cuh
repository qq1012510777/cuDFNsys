///////////////////////////////////////////////////////////////////
// NAME:              IntersectionTwo2DSegs.cuh
//
// PURPOSE:           Identify intersection between two 2D segments
//
// FUNCTIONS/OBJECTS: IntersectionTwo2DSegs
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../../GlobalDef/GlobalDef.cuh"
#include "If2DPntLiesOnCollinearSeg.cuh"
#include "IntersectionTwoCollinearSegs.cuh"
#include "OrientationThree2DPnts.cuh"

namespace cuDFNsys
{
__device__ __host__ bool IntersectionTwo2DSegs(float2 *Seg_1,
                                               float2 *Seg_2,
                                               float2 *intersection,
                                               int *sign, // 1, pnt; 2, seg; 3, none
                                               float _TOL_)
{
    bool state_ = false;

    float2 p1, q1, p2, q2;
    p1 = Seg_1[0];
    q1 = Seg_1[1];
    p2 = Seg_2[0];
    q2 = Seg_2[1];

    int o1 = cuDFNsys::OrientationThree2DPnts(p1, q1, p2);
    int o2 = cuDFNsys::OrientationThree2DPnts(p1, q1, q2);
    int o3 = cuDFNsys::OrientationThree2DPnts(p2, q2, p1);
    int o4 = cuDFNsys::OrientationThree2DPnts(p2, q2, q1);

    if (o1 == 0 && o2 == 0 && o3 == 0 && o4 == 0) // two line segments may overlap
    {
        int sign_2;
        float2 Seg_a[2] = {Seg_1[0], Seg_1[1]};
        float2 Seg_b[2] = {Seg_2[0], Seg_2[1]};
        float2 intersection_a[2];

        state_ = cuDFNsys::IntersectionTwoCollinearSegs(Seg_a, Seg_b, intersection_a, &sign_2, _TOL_);

        *sign = sign_2;
        intersection[0] = intersection_a[0];
        intersection[1] = intersection_a[1];
        return state_;
    };

    // General case
    if (o1 != o2 && o3 != o4)
        state_ = true;

    if (o1 == 0 && cuDFNsys::If2DPntLiesOnCollinearSeg(p1, p2, q1)) // Special Cases, p1, q1 and p2 are collinear and p2 lies on segment p1q1
    {
        *sign = 1;
        intersection[0] = p2;
        return true;
    }
    else if (o2 == 0 && cuDFNsys::If2DPntLiesOnCollinearSeg(p1, q2, q1)) // p1, q1 and q2 are collinear and q2 lies on segment p1q1
    {
        *sign = 1;
        intersection[0] = q2;
        return true;
    }
    else if (o3 == 0 && cuDFNsys::If2DPntLiesOnCollinearSeg(p2, p1, q2)) // p2, q2 and p1 are collinear and p1 lies on segment p2q2
    {
        *sign = 1;
        intersection[0] = p1;
        return true;
    }
    else if (o4 == 0 && cuDFNsys::If2DPntLiesOnCollinearSeg(p2, q1, q2)) // p2, q2 and q1 are collinear and q1 lies on segment p2q2
    {
        *sign = 1;
        intersection[0] = q1;
        return true;
    }

    if (state_ == false)
        return state_;

    // two line segments intersects in a general way

    // first, find intersection between two infinite lines

    // Line AB represented as a1x + b1y = c1
    float a1 = Seg_1[1].y - Seg_1[0].y;
    float b1 = Seg_1[0].x - Seg_1[1].x;
    float c1 = a1 * (Seg_1[0].x) + b1 * (Seg_1[0].y);

    // Line CD represented as a2x + b2y = c2
    float a2 = Seg_2[1].y - Seg_2[0].y;
    float b2 = Seg_2[0].x - Seg_2[1].x;
    float c2 = a2 * (Seg_2[0].x) + b2 * (Seg_2[0].y);

    float determinant = a1 * b2 - a2 * b1;

    intersection[0].x = (b2 * c1 - b1 * c2) / determinant;
    intersection[0].y = (a1 * c2 - a2 * c1) / determinant;

    *sign = 1;

    return true;
};
}; // namespace cuDFNsys