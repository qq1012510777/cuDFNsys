///////////////////////////////////////////////////////////////////
// NAME:              IfParticlePositionInNeighboringElement.cuh
//
// PURPOSE:           If the particle position is in a neighboring element
//
// FUNCTIONS/OBJECTS: IfParticlePositionInNeighboringElement
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
__device__ __host__ bool IfParticlePositionInNeighboringElement(float2 Position_p,
                                                                uint &EleID,
                                                                uint *EleToFracID_ptr,
                                                                uint FracID,
                                                                cuDFNsys::NeighborEle NE,
                                                                cuDFNsys::EleCoor *Coordinate2D_Vec_dev_ptr,
                                                                float _TOL_p)
{

    //---------- if it is inside a grid?
    //---------- if it is inside a grid?
    //---------- if it is inside a grid?
    for (uint j = 0; j < NE.NumNeighborEle; ++j)
    {
        uint eleID2 = NE.EleID[j];

        // for (uint p = 0; p < 10; ++p)
        //     if (VisitedEleID[p] != -1)
        //         if (eleID2 == (uint)VisitedEleID[p])
        //             continue;
        if (EleToFracID_ptr[eleID2 - 1] == FracID)
        {
            float2 Vertex_Triangle2[3];
            Vertex_Triangle2[0] = make_float2(Coordinate2D_Vec_dev_ptr[eleID2 - 1].x[0], Coordinate2D_Vec_dev_ptr[eleID2 - 1].y[0]);
            Vertex_Triangle2[1] = make_float2(Coordinate2D_Vec_dev_ptr[eleID2 - 1].x[1], Coordinate2D_Vec_dev_ptr[eleID2 - 1].y[1]);
            Vertex_Triangle2[2] = make_float2(Coordinate2D_Vec_dev_ptr[eleID2 - 1].x[2], Coordinate2D_Vec_dev_ptr[eleID2 - 1].y[2]);

            bool IfInGrid = cuDFNsys::IfPntInside2DConvexPoly(Position_p,
                                                              Vertex_Triangle2,
                                                              3);

            if (IfInGrid)
            {
                EleID = eleID2;
                return true;
            }
        }
    }

    //---------- if it is on one of the bounds of a grid
    //---------- if it is on one of the bounds of a grid
    //---------- if it is on one of the bounds of a grid
    for (uint j = 0; j < NE.NumNeighborEle; ++j)
    {
        uint eleID2 = NE.EleID[j];

        // for (uint p = 0; p < 10; ++p)
        //     if (VisitedEleID[p] != -1)
        //         if (eleID2 == (uint)VisitedEleID[p])
        //             continue;
        if (EleToFracID_ptr[eleID2 - 1] == FracID)
        {
            float2 Vertex_Triangle2[3];
            Vertex_Triangle2[0] = make_float2(Coordinate2D_Vec_dev_ptr[eleID2 - 1].x[0], Coordinate2D_Vec_dev_ptr[eleID2 - 1].y[0]);
            Vertex_Triangle2[1] = make_float2(Coordinate2D_Vec_dev_ptr[eleID2 - 1].x[1], Coordinate2D_Vec_dev_ptr[eleID2 - 1].y[1]);
            Vertex_Triangle2[2] = make_float2(Coordinate2D_Vec_dev_ptr[eleID2 - 1].x[2], Coordinate2D_Vec_dev_ptr[eleID2 - 1].y[2]);

            bool IfNearVertex = cuDFNsys::IfParticlePositionNearOneVertexOfElement(Position_p,
                                                                                   Vertex_Triangle2,
                                                                                   _TOL_p);

            int edgenolocal = -1;
            bool IfOnGridBound = cuDFNsys::IfParticleOnBoundOfElement(Position_p,
                                                                      Vertex_Triangle2,
                                                                      edgenolocal,
                                                                      _TOL_p);

            if (IfNearVertex == true ||
                IfOnGridBound == true)
            {
                EleID = eleID2;
                return true;
            }
        }
    }

    return false;
};
}; // namespace cuDFNsys