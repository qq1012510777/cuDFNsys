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
#include "../DataTypeSelector/DataTypeSelector.cuh"
#include "../GlobalDef/GlobalDef.cuh"
#include "../MatrixManipulation/MatrixManipulation.cuh"
#include "Fracture.cuh"

namespace cuDFNsys
{
// Truncate a fracture in a DFN
template <typename T>
__device__ __host__ bool TruncateFracture(cuDFNsys::Fracture<T> *verts,
                                          cuDFNsys::Vector1<T> L,
                                          int plane,
                                          int dir);
}; // namespace cuDFNsys