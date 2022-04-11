///////////////////////////////////////////////////////////////////
// NAME:              IdentifyFracPairSphericalDetection.cuh
//
// PURPOSE:           Identify fracture pair where the two circumscribed
//                    spheres intersect
//
// FUNCTIONS/OBJECTS: IdentifyFracPairSphericalDetection
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../GlobalDef/GlobalDef.cuh"
#include "Fracture.cuh"

namespace cuDFNsys
{
__global__ void IdentifyFracPairSphericalDetection(cuDFNsys::Fracture *verts,
                                                   int3 *Frac_pairs,
                                                   int InitialPairNO,
                                                   int count);
}; // namespace cuDFNsys