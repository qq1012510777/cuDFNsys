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
#include "../Geometry/Geometry.cuh"
#include "../Mesh/EleCoor.cuh"
#include "EdgeToEle.cuh"
#include "IfParticlePositionNearOneVertexOfElement.cuh"
#include "NeighborEle.cuh"

namespace cuDFNsys
{
template <typename T>
__device__ __host__ bool IfParticlePositionInNeighboringElement(cuDFNsys::Vector2<T> Position_p,
                                                                uint &EleID,
                                                                uint *EleToFracID_ptr,
                                                                uint FracID,
                                                                cuDFNsys::NeighborEle NE,
                                                                cuDFNsys::EleCoor<T> *Coordinate2D_Vec_dev_ptr,
                                                                T _TOL_p)
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
            cuDFNsys::Vector2<T> Vertex_Triangle2[3];
            Vertex_Triangle2[0] = cuDFNsys::MakeVector2(Coordinate2D_Vec_dev_ptr[eleID2 - 1].x[0], Coordinate2D_Vec_dev_ptr[eleID2 - 1].y[0]);
            Vertex_Triangle2[1] = cuDFNsys::MakeVector2(Coordinate2D_Vec_dev_ptr[eleID2 - 1].x[1], Coordinate2D_Vec_dev_ptr[eleID2 - 1].y[1]);
            Vertex_Triangle2[2] = cuDFNsys::MakeVector2(Coordinate2D_Vec_dev_ptr[eleID2 - 1].x[2], Coordinate2D_Vec_dev_ptr[eleID2 - 1].y[2]);

            bool IfInGrid = cuDFNsys::IfPntInside2DConvexPoly<T>(Position_p,
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
            cuDFNsys::Vector2<T> Vertex_Triangle2[3];
            Vertex_Triangle2[0] = cuDFNsys::MakeVector2(Coordinate2D_Vec_dev_ptr[eleID2 - 1].x[0], Coordinate2D_Vec_dev_ptr[eleID2 - 1].y[0]);
            Vertex_Triangle2[1] = cuDFNsys::MakeVector2(Coordinate2D_Vec_dev_ptr[eleID2 - 1].x[1], Coordinate2D_Vec_dev_ptr[eleID2 - 1].y[1]);
            Vertex_Triangle2[2] = cuDFNsys::MakeVector2(Coordinate2D_Vec_dev_ptr[eleID2 - 1].x[2], Coordinate2D_Vec_dev_ptr[eleID2 - 1].y[2]);

            bool IfNearVertex = cuDFNsys::IfParticlePositionNearOneVertexOfElement<T>(Position_p,
                                                                                      Vertex_Triangle2,
                                                                                      _TOL_p);

            int edgenolocal = -1;
            bool IfOnGridBound = cuDFNsys::IfParticleOnBoundOfElement<T>(Position_p,
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