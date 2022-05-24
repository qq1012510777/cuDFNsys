///////////////////////////////////////////////////////////////////
// NAME:              IfParticleOnBoundOfElement.cuh
//
// PURPOSE:           If particle lies on bound of element
//
// FUNCTIONS/OBJECTS: IfParticleOnBoundOfElement
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
#include "IfParticlePositionNearOneVertexOfElement.cuh"
#include "Particle.cuh"

namespace cuDFNsys
{
__host__ __device__ bool IfParticleOnBoundOfElement(float2 PositionP,
                                                    float2 GridVertex[3],
                                                    int &EdgeNOLocal,
                                                    float _TOL_)
{

    for (uint i = 0; i < 3; ++i)
    {
        int o1 = cuDFNsys::OrientationThree2DPnts(GridVertex[i],
                                                  GridVertex[(i + 1) % 3],
                                                  PositionP, _TOL_);
        if (o1 == 0 && cuDFNsys::If2DPntLiesOnCollinearSeg(GridVertex[i],
                                                           PositionP,
                                                           GridVertex[(i + 1) % 3]))
        {
            EdgeNOLocal = (int)i;
            return true;
        };
    };

    EdgeNOLocal = -1;
    return false;
};
}; // namespace cuDFNsys