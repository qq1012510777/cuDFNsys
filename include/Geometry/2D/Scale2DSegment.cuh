///////////////////////////////////////////////////////////////////
// NAME:              Scale2DSegment.cuh
//
// PURPOSE:           scale a 2D line segment along the center
//
// FUNCTIONS/OBJECTS: Scale2DSegment
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../../GlobalDef/GlobalDef.cuh"
#include "../../DataTypeSelector/DataTypeSelector.cuh"
namespace cuDFNsys
{
template <typename T>
__host__ __device__ void Scale2DSegment(cuDFNsys::Vector2<T> *Segment, T scaleF);
}; // namespace cuDFNsys