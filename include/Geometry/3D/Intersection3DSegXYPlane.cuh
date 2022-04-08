///////////////////////////////////////////////////////////////////
// NAME:              Intersection3DSegXYPlane.cuh
//
// PURPOSE:           Identify intersection between
//                    a 3D segment and the XY plane
//
// FUNCTIONS/OBJECTS: Intersection3DSegXYPlane
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../../GlobalDef/GlobalDef.cuh"

namespace cuDFNsys
{
// intersection between a 3D segment and the XY plane
__device__ __host__ bool Intersection3DSegXYPlane(float3 *Seg,
                                                  float3 *Intersec_PNT,
                                                  int *sign_,
                                                  float _TOL_); // sign: 1: pnt; 2: seg; -1: none;
};                                                              // namespace cuDFNsys