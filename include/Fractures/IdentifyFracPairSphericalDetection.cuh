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
#include "../DataTypeSelector/DataTypeSelector.cuh"
#include "../GlobalDef/GlobalDef.cuh"
#include "../MatrixManipulation/MatrixManipulation.cuh"
#include "Fracture.cuh"

namespace cuDFNsys
{
template <typename T>
__global__ void IdentifyFracPairSphericalDetection(cuDFNsys::Fracture<T> *verts,
                                                   int3 *Frac_pairs,
                                                   int InitialPairNO,
                                                   int count);
}; // namespace cuDFNsys