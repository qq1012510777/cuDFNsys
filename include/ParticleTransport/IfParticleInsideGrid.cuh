///////////////////////////////////////////////////////////////////
// NAME:              IfParticleInsideGrid.cuh
//
// PURPOSE:           If the next position of particle is still inside the
//                    grid
//
// FUNCTIONS/OBJECTS: IfParticleInsideGrid
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
__device__ __host__ bool IfParticleInsideGrid(float2 Position_p, float2 Vertex_Triangle[3])
{
    bool InGrid = cuDFNsys::IfPntInside2DConvexPoly(Position_p,
                                                    Vertex_Triangle, 3);

    bool OnGridBound = cuDFNsys::IfPntLiesOnBound2DConvexPoly(Position_p,
                                                              Vertex_Triangle,
                                                              3,
                                                              _TOL_ParticleOnGridBound);

    if (InGrid == true || OnGridBound == true)
        return true;

    return false;
};
}; // namespace cuDFNsys