///////////////////////////////////////////////////////////////////
// NAME:              Triangle2DArea.cuh
//
// PURPOSE:           Get 2D triangle area
//
// FUNCTIONS/OBJECTS: Triangle2DArea
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../../DataTypeSelector/DataTypeSelector.cuh"
#include "../../GlobalDef/GlobalDef.cuh"
#include "../../MatrixManipulation/MatrixManipulation.cuh"

namespace cuDFNsys
{
template <typename T>
__device__ __host__ T Triangle2DArea(cuDFNsys::Vector2<T> Pnt1,
                                     cuDFNsys::Vector2<T> Pnt2,
                                     cuDFNsys::Vector2<T> Pnt3);
}; // namespace cuDFNsys