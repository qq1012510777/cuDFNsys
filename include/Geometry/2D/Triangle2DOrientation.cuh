///////////////////////////////////////////////////////////////////
// NAME:              Triangle2DOrientation.cuh
//
// PURPOSE:           return orientation of a 2D triangle
//                    true: clockwise
//                    false: counterclockwise
//
// FUNCTIONS/OBJECTS: Triangle2DOrientation
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../../DataTypeSelector/DataTypeSelector.cuh"
#include "../../GlobalDef/GlobalDef.cuh"

namespace cuDFNsys
{
template <typename T>
__device__ __host__ bool Triangle2DOrientation(cuDFNsys::Vector2<T> a,
                                               cuDFNsys::Vector2<T> b,
                                               cuDFNsys::Vector2<T> c);
}; // namespace cuDFNsys