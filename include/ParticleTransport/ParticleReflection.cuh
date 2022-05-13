///////////////////////////////////////////////////////////////////
// NAME:              ParticleReflection.cuh
//
// PURPOSE:           Reflect a particle if the particle crosses a non-flux
//                    edge
//
// FUNCTIONS/OBJECTS: ParticleReflection
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
__device__ __host__ float2 ParticleReflection(float2 P, float2 A, float2 B)
{
    // Performing translation and shifting origin at A
    float2 Pt = make_float2(P.x - A.x, P.y - A.y);
    float2 Bt = make_float2(B.x - A.x, B.y - A.y);

    // Performing rotation in clockwise direction
    // BtAt becomes the X-Axis in the new coordinate system
    float2 Pr = make_float2((Pt.x * Bt.x + Pt.y * Bt.y) / (Bt.x * Bt.x + Bt.y * Bt.y),
                            -1.0f * (Pt.y * Bt.x - Pt.x * Bt.y) / (Bt.x * Bt.x + Bt.y * Bt.y));

    //printf("Pr: %f %f\n", Pr.x, Pr.y);
    // Reflection of Pr about the new X-Axis
    // Followed by restoring from rotation
    // Followed by restoring from translation

    float2 Ps;
    Ps.x = Pr.x * Bt.x - Pr.y * Bt.y + A.x;
    Ps.y = Pr.x * Bt.y + Pr.y * Bt.x + A.y;

    return Ps;
};
}; // namespace cuDFNsys