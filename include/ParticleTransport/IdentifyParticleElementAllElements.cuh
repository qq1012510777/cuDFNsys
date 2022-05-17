///////////////////////////////////////////////////////////////////
// NAME:              IdentifyParticleElementAllElements.cuh
//
// PURPOSE:           Identify which element does the particle lies in
//
// FUNCTIONS/OBJECTS: IdentifyParticleElementAllElements
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../Geometry/2D/DistancePnt2DSeg.cuh"
#include "../Geometry/2D/IfPntInside2DConvexPoly.cuh"
#include "../Geometry/2D/IfPntLiesOnBound2DConvexPoly.cuh"
#include "../Geometry/2D/IntersectionTwo2DSegs.cuh"
#include "../Geometry/2D/Scale2DSegment.cuh"
#include "../MHFEM/MHFEM.cuh"
#include "../MHFEM/ReconstructVelocityGrid.cuh"
#include "EdgeToEle.cuh"
#include "IfParticleInsideGrid.cuh"
#include "Particle.cuh"

namespace cuDFNsys
{
__host__ __device__ bool IdentifyParticleElementAllElements(float2 Position_p1,
                                                            float2 Position_p2,
                                                            uint &EleID,
                                                            int numElements,
                                                            uint *EleToFracID_ptr,
                                                            uint FracID,
                                                            cuDFNsys::NeighborEle NE,
                                                            cuDFNsys::EleCoor *Coordinate2D_Vec_dev_ptr,
                                                            int VisitedEleID[10])
{
    int eleID_end1LiesIn = -1;

    for (uint j = 0; j < NE.NumNeighborEle; ++j)
    {
        uint eleID2 = NE.EleID[j];

        for (uint p = 0; p < 10; ++p)
            if (VisitedEleID[p] != -1)
                if (eleID2 == (uint)VisitedEleID[p])
                    continue;

        float2 Vertex_Triangle2[3];
        Vertex_Triangle2[0] = make_float2(Coordinate2D_Vec_dev_ptr[eleID2 - 1].x[0], Coordinate2D_Vec_dev_ptr[eleID2 - 1].y[0]);
        Vertex_Triangle2[1] = make_float2(Coordinate2D_Vec_dev_ptr[eleID2 - 1].x[1], Coordinate2D_Vec_dev_ptr[eleID2 - 1].y[1]);
        Vertex_Triangle2[2] = make_float2(Coordinate2D_Vec_dev_ptr[eleID2 - 1].x[2], Coordinate2D_Vec_dev_ptr[eleID2 - 1].y[2]);

        bool InGrid2 = cuDFNsys::IfParticleInsideGrid(Position_p1,
                                                      Vertex_Triangle2);
        bool InGrid3 = cuDFNsys::IfParticleInsideGrid(Position_p2,
                                                      Vertex_Triangle2);

        bool InGrid4 = cuDFNsys::IfPntInside2DConvexPoly(Position_p1,
                                                         Vertex_Triangle2, 3);
        if (InGrid2 == true && InGrid3 == true)
        {
            // we do not need to change grid ID in this case
            EleID = eleID2;
            //printf("In neighboring elements_\n");
            return true;
        };

        if (InGrid4 == true)
        {
            eleID_end1LiesIn = (int)eleID2;
        };
    }

    for (int k = 0; k < numElements; ++k)
    {

        if (EleToFracID_ptr[k] == FracID)
        {
            for (uint p = 0; p < 10; ++p)
                if (VisitedEleID[p] != -1)
                    if ((k + 1) == (uint)VisitedEleID[p])
                        continue;

            float2 Vertex_Triangle3[3];
            Vertex_Triangle3[0] = make_float2(Coordinate2D_Vec_dev_ptr[k].x[0], Coordinate2D_Vec_dev_ptr[k].y[0]);
            Vertex_Triangle3[1] = make_float2(Coordinate2D_Vec_dev_ptr[k].x[1], Coordinate2D_Vec_dev_ptr[k].y[1]);
            Vertex_Triangle3[2] = make_float2(Coordinate2D_Vec_dev_ptr[k].x[2], Coordinate2D_Vec_dev_ptr[k].y[2]);

            bool InGrid2 = cuDFNsys::IfParticleInsideGrid(Position_p1,
                                                          Vertex_Triangle3);
            bool InGrid3 = cuDFNsys::IfParticleInsideGrid(Position_p2,
                                                          Vertex_Triangle3);

            bool InGrid4 = cuDFNsys::IfPntInside2DConvexPoly(Position_p1,
                                                             Vertex_Triangle3, 3);
            if (InGrid2 == true && InGrid3 == true)
            {
                EleID = k + 1;
                return true;
            }

            if (InGrid4 == true)
            {
                eleID_end1LiesIn = (int)(k + 1);
            }
        }
    }

    if (eleID_end1LiesIn != -1)
    {
        for (uint p = 0; p < 10; ++p)
            if (VisitedEleID[p] != -1)
                if (eleID_end1LiesIn == (uint)VisitedEleID[p])
                    return false;

        EleID = (uint)eleID_end1LiesIn;

        return true;
    }

    return false;
}
}; // namespace cuDFNsys