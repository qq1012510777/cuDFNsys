///////////////////////////////////////////////////////////////////
// NAME:              IdentifyIntersectionKernel.cuh
//
// PURPOSE:           Identify intersection using GPU
//
// FUNCTIONS/OBJECTS: IdentifyIntersectionKernel
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../Geometry/Geometry.cuh"
#include "../GlobalDef/GlobalDef.cuh"
#include "../MatrixManipulation/MatrixManipulation.cuh"
#include "Fracture.cuh"
#include "Intersection.cuh"
using namespace std;

namespace cuDFNsys
{
template <typename T>
__global__ void IdentifyIntersectionKernel(cuDFNsys::Fracture<T> *verts,
                                           int count,
                                           cuDFNsys::Intersection<T> *Int_sec,
                                           bool if_trucncated);
}; // namespace cuDFNsys