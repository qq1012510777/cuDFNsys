///////////////////////////////////////////////////////////////////
// NAME:              Roate2DPositionTo3D.cuh
//
// PURPOSE:           Rotate a 2D particle position to 3D
//
// FUNCTIONS/OBJECTS: Roate2DPositionTo3D
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
__host__ __device__ float3 Roate2DPositionTo3D(float2 PositionP,
                                               cuDFNsys::Fracture OneFrac)
{
    float3 P_3D = make_float3(PositionP.x, PositionP.y, 0.0f);

    float RK_2[3][3];
    OneFrac.RoationMatrix(RK_2, 23);

    P_3D = cuDFNsys::ProductSquare3Float3(RK_2, P_3D);
    P_3D.x += OneFrac.Center.x;
    P_3D.y += OneFrac.Center.y;
    P_3D.z += OneFrac.Center.z;

    return P_3D;
}
}; // namespace cuDFNsys