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
#include "../../GlobalDef/GlobalDef.cuh"
#include "DistancePnt2DSeg.cuh"

namespace cuDFNsys
{
//if a point lies on 2D bound
__device__ __host__ bool IfPntLiesOnBound2DConvexPoly(float2 pnt,
                                                      float2 *verts,
                                                      int N,
                                                      float _tol_);
// if a point lies on 2D bound, also return the edge NO
__device__ __host__ bool IfPntLiesOnBound2DConvexPolyReturnEdgeNO(float2 pnt,
                                                                  float2 *verts,
                                                                  int N,
                                                                  float _tol_,
                                                                  int *edgeNO);
}; // namespace cuDFNsys