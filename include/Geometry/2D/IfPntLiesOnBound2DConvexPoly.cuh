///////////////////////////////////////////////////////////////////
// NAME:              IfPntLiesOnBound2DConvexPoly.cuh
//
// PURPOSE:           if a point lies on bound of 2D convex polygon 2D
//
// FUNCTIONS/OBJECTS: IfPntLiesOnBound2DConvexPoly
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../../DataTypeSelector/DataTypeSelector.cuh"
#include "../../GlobalDef/GlobalDef.cuh"
#include "../../MatrixManipulation/MatrixManipulation.cuh"
#include "DistancePnt2DSeg.cuh"

namespace cuDFNsys
{
//if a point lies on 2D bound
template <typename T>
__device__ __host__ bool IfPntLiesOnBound2DConvexPoly(cuDFNsys::Vector2<T> pnt,
                                                      cuDFNsys::Vector2<T> *verts,
                                                      int N,
                                                      T _tol_);

// if a point lies on 2D bound, also return the edge NO
template <typename T>
__device__ __host__ bool IfPntLiesOnBound2DConvexPolyReturnEdgeNO(cuDFNsys::Vector2<T> pnt,
                                                                  cuDFNsys::Vector2<T> *verts,
                                                                  int N,
                                                                  T _tol_,
                                                                  int *edgeNO);
}; // namespace cuDFNsys