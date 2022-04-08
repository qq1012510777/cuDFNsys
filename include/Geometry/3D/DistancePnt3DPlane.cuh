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

namespace cuDFNsys
{
__device__ __host__ float DistancePnt3DPlane(float3 Plane[3], float3 pnt);
}; // namespace cuDFNsys