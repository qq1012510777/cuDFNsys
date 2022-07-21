///////////////////////////////////////////////////////////////////
// NAME:              DistancePnt3DPlane.cuh
//
// PURPOSE:           Distance between a 3D point and a 3D plane
//
// FUNCTIONS/OBJECTS: DistancePnt3DPlane
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../../GlobalDef/GlobalDef.cuh"
#include "../../DataTypeSelector/DataTypeSelector.cuh"

namespace cuDFNsys
{
template <typename T>
__device__ __host__ T DistancePnt3DPlane(cuDFNsys::Vector3<T> Plane[3], cuDFNsys::Vector3<T> pnt);
}; // namespace cuDFNsys