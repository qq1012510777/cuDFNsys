#include "Geometry/2D/Scale2DSegment.cuh"

// ====================================================
// NAME:        Scale2DSegment
// DESCRIPTION: scale a 2D line segment along the center
// AUTHOR:      Tingchang YIN
// DATE:        09/05/2022
// ====================================================
__host__ __device__ void cuDFNsys::Scale2DSegment(float2 *Segment, float scaleF)
{
    float2 center_;
    center_.x = 0.5f * (Segment[0].x + Segment[1].x);
    center_.y = 0.5f * (Segment[0].y + Segment[1].y);

    Segment[0].x = center_.x + (Segment[0].x - center_.x) * scaleF;
    Segment[0].y = center_.y + (Segment[0].y - center_.y) * scaleF;

    Segment[1].x = center_.x + (Segment[1].x - center_.x) * scaleF;
    Segment[1].y = center_.y + (Segment[1].y - center_.y) * scaleF;
}; // Scale2DSegment