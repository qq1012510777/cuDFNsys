///////////////////////////////////////////////////////////////////
// NAME:              Intersection3DPolyXYPlane.cuh
//
// PURPOSE:           Identify intersection between
//                    a 3D Polygon and the XY plane
//
// FUNCTIONS/OBJECTS: Intersection3DPolyXYPlane
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../../GlobalDef/GlobalDef.cuh"
#include "Intersection3DSegXYPlane.cuh"

namespace cuDFNsys
{
// Identify intersection between a 3D Polygon and the XY plane
__device__ __host__ bool Intersection3DPolyXYPlane(float3 *Poly,
                                                   int NUM_vert,
                                                   float3 *Intersection,
                                                   float _TOL_);
}; // namespace cuDFNsys