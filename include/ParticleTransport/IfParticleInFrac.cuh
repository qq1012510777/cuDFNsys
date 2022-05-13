///////////////////////////////////////////////////////////////////
// NAME:              IfParticleInFrac.cuh
//
// PURPOSE:           If the next position of particle is still inside the
//                    fracture
//
// FUNCTIONS/OBJECTS: IfParticleInFrac
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
__device__ __host__ bool IfParticleInFrac(float2 Position_p, float2 *vertex_frac2D, int NumVertexes)
{
    bool Infrac = cuDFNsys::IfPntInside2DConvexPoly(Position_p,
                                                    vertex_frac2D, NumVertexes);

    bool OnfracBound = cuDFNsys::IfPntLiesOnBound2DConvexPoly(Position_p,
                                                              vertex_frac2D,
                                                              NumVertexes, _TOL_ParticleOnGridBound);
    if (Infrac == false && OnfracBound == false)
        return false;
    return true;
};
}; // namespace cuDFNsys