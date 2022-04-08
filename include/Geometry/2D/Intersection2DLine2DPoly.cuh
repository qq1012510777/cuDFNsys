///////////////////////////////////////////////////////////////////
// NAME:              Intersection2DLine2DPoly.cuh
//
// PURPOSE:           Intersection between
//                    2D Line (infinite) and 2D Polygon
//
// FUNCTIONS/OBJECTS: Intersection2DLine2DPoly
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////

#pragma once
#include "../../GlobalDef/GlobalDef.cuh"
#include "Intersection2DLine2DSeg.cuh"

namespace cuDFNsys
{
// Intersection between 2D Line (infinite) and 2D Polygon
__device__ __host__ bool Intersection2DLine2DPoly(float2 *Poly2D,
                                                  int NUM_verts,
                                                  float2 *Line,
                                                  float2 *intersection_k,
                                                float _TOL_);
}; // namespace cuDFNsys