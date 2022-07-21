///////////////////////////////////////////////////////////////////
// NAME:              If2DPntLiesOnCollinearSeg.cuh
//
// PURPOSE:           check if a 2D point lies on a collinear segment
//                    point q, line segment: p-r
//
// FUNCTIONS/OBJECTS: If2DPntLiesOnCollinearSeg
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../../DataTypeSelector/DataTypeSelector.cuh"
#include "../../GlobalDef/GlobalDef.cuh"

namespace cuDFNsys
{
template <typename T>
__device__ __host__ bool If2DPntLiesOnCollinearSeg(cuDFNsys::Vector2<T> p, cuDFNsys::Vector2<T> q, cuDFNsys::Vector2<T> r);
}; // namespace cuDFNsys