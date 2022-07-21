///////////////////////////////////////////////////////////////////
// NAME:              ReconstructVelocityGrid.cuh
//
// PURPOSE:           Reconstruct the velocity field in a grid
//
// FUNCTIONS/OBJECTS: ReconstructVelocityGrid
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../DataTypeSelector/DataTypeSelector.cuh"
#include "../Geometry/Geometry.cuh"

namespace cuDFNsys
{
template <typename T>
__host__ __device__ cuDFNsys::Vector2<T> ReconstructVelocityGrid(cuDFNsys::Vector2<T> Point_,
                                                                 cuDFNsys::Vector2<T> Vertex[3],
                                                                 cuDFNsys::Vector3<T> VelocityEdgeNormal);
}; // namespace cuDFNsys