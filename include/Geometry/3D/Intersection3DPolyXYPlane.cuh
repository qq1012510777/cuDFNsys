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
#include "../../DataTypeSelector/DataTypeSelector.cuh"
#include "../../GlobalDef/GlobalDef.cuh"
#include "../../MatrixManipulation/MatrixManipulation.cuh"
#include "Intersection3DSegXYPlane.cuh"

namespace cuDFNsys
{
// Identify intersection between a 3D Polygon and the XY plane
template <typename T>
__device__ __host__ bool Intersection3DPolyXYPlane(cuDFNsys::Vector3<T> *Poly,
                                                   int NUM_vert,
                                                   cuDFNsys::Vector3<T> *Intersection,
                                                   T _TOL_);
}; // namespace cuDFNsys