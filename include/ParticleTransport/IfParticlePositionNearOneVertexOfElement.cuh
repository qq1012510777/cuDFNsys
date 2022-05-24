///////////////////////////////////////////////////////////////////
// NAME:              IfParticlePositionNearOneVertexOfElement.cuh
//
// PURPOSE:           If the particle position is very close to one of the vertexes of grid
//
// FUNCTIONS/OBJECTS: IfParticlePositionNearOneVertexOfElement
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../Geometry/2D/DistancePnt2DSeg.cuh"
#include "../Geometry/2D/IfPntInside2DConvexPoly.cuh"
#include "../Geometry/2D/IfPntLiesOnBound2DConvexPoly.cuh"
#include "../Geometry/2D/IntersectionTwo2DSegs.cuh"
#include "../MHFEM/MHFEM.cuh"
#include "../MHFEM/ReconstructVelocityGrid.cuh"
#include "EdgeToEle.cuh"
#include "Particle.cuh"

namespace cuDFNsys
{
__device__ __host__ bool IfParticlePositionNearOneVertexOfElement(float2 Position_p, float2 Tri[3], float _TOL_p)
{
    for (uint i = 0; i < 3; ++i)
    {
        float2 V;
        V.x = Position_p.x - Tri[i].x;
        V.y = Position_p.y - Tri[i].y;

        float norm_ = sqrt(V.x * V.x + V.y * V.y);

        if (norm_ < _TOL_p)
            return true;
    }
    return false;
};
}; // namespace cuDFNsys