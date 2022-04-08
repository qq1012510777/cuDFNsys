///////////////////////////////////////////////////////////////////
// NAME:              TruncateFracture.cuh
//
// PURPOSE:           Truncate a fracture in a DFN,
//                    function as 3D polygon clipping
//                    by a box
//
// FUNCTIONS/OBJECTS: Fractures
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../GlobalDef/GlobalDef.cuh"
#include "Fracture.cuh"

namespace cuDFNsys
{
// Truncate a fracture in a DFN
__device__ __host__ bool TruncateFracture(cuDFNsys::Fracture *verts,
                                          float L,
                                          int plane,
                                          int dir);
}; // namespace cuDFNsys