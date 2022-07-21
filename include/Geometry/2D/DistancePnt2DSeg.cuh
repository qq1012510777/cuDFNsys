///////////////////////////////////////////////////////////////////
// NAME:              DistancePnt2DSeg.cuh
//
// PURPOSE:           get distance between 2D point and 2D segment
//
// FUNCTIONS/OBJECTS: DistancePnt2DSeg
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
__device__ __host__ T DistancePnt2DSeg(cuDFNsys::Vector2<T> pnt, cuDFNsys::Vector2<T> *verts);
}; // namespace cuDFNsys