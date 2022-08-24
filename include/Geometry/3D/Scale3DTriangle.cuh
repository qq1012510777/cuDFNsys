///////////////////////////////////////////////////////////////////
// NAME:              Scale3DTriangle.cuh
//
// PURPOSE:           Scale a 3D triangle
//
// FUNCTIONS/OBJECTS: Scale3DTriangle
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
__device__ __host__ void Scale3DTriangle(cuDFNsys::Vector3<T> Triangle[3],
                                         T factor);
}; // namespace cuDFNsys