///////////////////////////////////////////////////////////////////
// NAME:              IntersectionTwoCollinear2DSegs.cuh
//
// PURPOSE:           Identify intersection between two collinear segments
//
// FUNCTIONS/OBJECTS: IntersectionTwoCollinear2DSegs
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../../DataTypeSelector/DataTypeSelector.cuh"
#include "../../GlobalDef/GlobalDef.cuh"
#include "../../MatrixManipulation/MatrixManipulation.cuh"
namespace cuDFNsys
{
template <typename T>
__device__ __host__ bool IntersectionTwoCollinear2DSegs(cuDFNsys::Vector2<T> *Seg_1,
                                                        cuDFNsys::Vector2<T> *Seg_2,
                                                        cuDFNsys::Vector2<T> *intersection,
                                                        int *sign,
                                                        T _TOL_);
}; // namespace cuDFNsys