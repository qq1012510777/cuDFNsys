///////////////////////////////////////////////////////////////////
// NAME:              Triangle3DArea.cuh
//
// PURPOSE:           Get 3D triangle area
//
// FUNCTIONS/OBJECTS: Triangle3DArea
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
__device__ __host__ T Triangle3DArea(cuDFNsys::Vector3<T> Pnt1,
                                     cuDFNsys::Vector3<T> Pnt2,
                                     cuDFNsys::Vector3<T> Pnt3);
}; // namespace cuDFNsys