///////////////////////////////////////////////////////////////////
// NAME:              Intersection2DLine2DSeg.cuh
//
// PURPOSE:           Identify intersection between
//                    2D line (infinite) and 2D segment
//
// FUNCTIONS/OBJECTS: Intersection2DLine2DSeg
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////

#pragma once
#include "../../DataTypeSelector/DataTypeSelector.cuh"
#include "../../GlobalDef/GlobalDef.cuh"
#include "../../MatrixManipulation/MatrixManipulation.cuh"

namespace cuDFNsys
{
// Identify intersection between 2D line (infinite) and 2D segment
template <typename T>
__device__ __host__ bool Intersection2DLine2DSeg(cuDFNsys::Vector2<T> *Line,
                                                 cuDFNsys::Vector2<T> *Seg,
                                                 int *sign_, // 1, pnt; 2, seg; 3, none
                                                 cuDFNsys::Vector2<T> *intersection,
                                                 T _TOL_);
}; // namespace cuDFNsys
