///////////////////////////////////////////////////////////////////
// NAME:              OrientationThree2DPnts.cuh
//
// PURPOSE:           return orientation of three 2D points
//                    0: collinear
//                    1: clockwise
//                    2: counterclockwise
//
// FUNCTIONS/OBJECTS: OrientationThree2DPnts
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../../DataTypeSelector/DataTypeSelector.cuh"
#include "../../GlobalDef/GlobalDef.cuh"

namespace cuDFNsys
{
template <typename T>
__device__ __host__ uint OrientationThree2DPnts(cuDFNsys::Vector2<T> p, cuDFNsys::Vector2<T> q, cuDFNsys::Vector2<T> r, T _tol_);
}; // namespace cuDFNsys