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
#include "../../DataTypeSelector/DataTypeSelector.cuh"
#include "../../GlobalDef/GlobalDef.cuh"
#include "../../MatrixManipulation/MatrixManipulation.cuh"

namespace cuDFNsys
{
// intersection between a 3D segment and the XY plane
template <typename T>
__device__ __host__ bool Intersection3DSegXYPlane(cuDFNsys::Vector3<T> *Seg,
                                                  cuDFNsys::Vector3<T> *Intersec_PNT,
                                                  int *sign_, // sign: 1: pnt; 2: seg; -1: none;
                                                  T _TOL_);
}; // namespace cuDFNsys