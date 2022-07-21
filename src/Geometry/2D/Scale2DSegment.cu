#include "Geometry/2D/Scale2DSegment.cuh"

// ====================================================
// NAME:        Scale2DSegment
// DESCRIPTION: scale a 2D line segment along the center
// AUTHOR:      Tingchang YIN
// DATE:        09/05/2022
// ====================================================
template <typename T>
__host__ __device__ void cuDFNsys::Scale2DSegment(cuDFNsys::Vector2<T> *Segment, T scaleF)
{
    cuDFNsys::Vector2<T> center_;
    center_.x = (T)0.5 * (Segment[0].x + Segment[1].x);
    center_.y = (T)0.5 * (Segment[0].y + Segment[1].y);

    Segment[0].x = center_.x + (Segment[0].x - center_.x) * scaleF;
    Segment[0].y = center_.y + (Segment[0].y - center_.y) * scaleF;

    Segment[1].x = center_.x + (Segment[1].x - center_.x) * scaleF;
    Segment[1].y = center_.y + (Segment[1].y - center_.y) * scaleF;
}; // Scale2DSegment
template __host__ __device__ void cuDFNsys::Scale2DSegment<double>(cuDFNsys::Vector2<double> *Segment, double scaleF);
template __host__ __device__ void cuDFNsys::Scale2DSegment<float>(cuDFNsys::Vector2<float> *Segment, float scaleF);