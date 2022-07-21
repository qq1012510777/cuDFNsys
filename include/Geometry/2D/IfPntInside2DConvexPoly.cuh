///////////////////////////////////////////////////////////////////
// NAME:              IfPntInside2DConvexPoly.cuh
//
// PURPOSE:           if a point is inside a 2D polygon
//
// FUNCTIONS/OBJECTS: IfPntInside2DConvexPoly
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../../GlobalDef/GlobalDef.cuh"
#include "../../MatrixManipulation/MatrixManipulation.cuh"

namespace cuDFNsys
{
template <typename T>
__device__ __host__ bool IfPntInside2DConvexPoly(cuDFNsys::Vector2<T> pnt,
                                                 cuDFNsys::Vector2<T> *verts,
                                                 int N);
}; // namespace cuDFNsys