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
#include "../../GlobalDef/GlobalDef.cuh"

namespace cuDFNsys
{
__device__ __host__ float Triangle3DArea(float3 Pnt1,
                                         float3 Pnt2,
                                         float3 Pnt3);
}; // namespace cuDFNsys