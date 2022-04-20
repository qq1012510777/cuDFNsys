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
#include "../Geometry/2D/Triangle2DArea.cuh"

namespace cuDFNsys
{
__host__ __device__ float2 ReconstructVelocityGrid(float2 Point_,
                                                   float2 Vertex[3],
                                                   float3 VelocityEdgeNormal);
}; // namespace cuDFNsys