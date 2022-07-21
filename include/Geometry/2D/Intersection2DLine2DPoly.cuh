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
#include "../../DataTypeSelector/DataTypeSelector.cuh"
#include "../../GlobalDef/GlobalDef.cuh"
#include "../../MatrixManipulation/MatrixManipulation.cuh"
#include "Intersection2DLine2DSeg.cuh"

namespace cuDFNsys
{
// Intersection between 2D Line (infinite) and 2D Polygon
template <typename T>
__device__ __host__ bool Intersection2DLine2DPoly(cuDFNsys::Vector2<T> *Poly2D,
                                                  int NUM_verts,
                                                  cuDFNsys::Vector2<T> *Line,
                                                  cuDFNsys::Vector2<T> *intersection_k,
                                                  T _TOL_);
}; // namespace cuDFNsys