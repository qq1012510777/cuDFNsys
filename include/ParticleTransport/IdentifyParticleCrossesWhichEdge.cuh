///////////////////////////////////////////////////////////////////
// NAME:              IdentifyParticleCrossesWhichEdge.cuh
//
// PURPOSE:           Identify which edge does the particle cross
//                      -1: no
//                       0 1 2
//
// FUNCTIONS/OBJECTS: IdentifyParticleCrossesWhichEdge
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
#include "Particle.cuh"

namespace cuDFNsys
{
__host__ __device__ float3 IdentifyParticleCrossesWhichEdge(float2 *P_Trajectory,
                                                            float2 *Vertex_Triangle,
                                                            float _TOL_,
                                                            float2 CrossedGlobalEdge[10][2],
                                                            int CountCrossedGlobalEdge,
                                                            uint stepNO,
                                                            uint particleNO)
{
    bool If_End1_on_Bound = false;
    int edge_onBound_end1 = -1;

    float3 Result;

    float2 p1, q1;
    p1 = P_Trajectory[0];
    q1 = P_Trajectory[1];

    for (uint i = 0; i < 3; ++i)
    {
        float2 p2, q2;
        p2 = Vertex_Triangle[i];
        q2 = Vertex_Triangle[(i + 1) % 3];

        // check if the edge has been crossed?
        // check if the edge has been crossed?
        // check if the edge has been crossed?
        if (CountCrossedGlobalEdge > 0)
        {
            bool IfEdgeChecked = false;
            for (int k = 0; k < CountCrossedGlobalEdge; ++k)
            {
                float2 a = CrossedGlobalEdge[k][0];
                float2 b = CrossedGlobalEdge[k][1];

                if (abs(a.x - p2.x) < 1e-7 && abs(a.y - p2.y) < 1e-7 &&
                    abs(b.x - q2.x) < 1e-7 && abs(b.y - q2.y) < 1e-7)
                {
                    IfEdgeChecked = true;
                    break;
                }

                if (abs(b.x - p2.x) < 1e-7 && abs(b.y - p2.y) < 1e-7 &&
                    abs(a.x - q2.x) < 1e-7 && abs(a.y - q2.y) < 1e-7)
                {
                    IfEdgeChecked = true;
                    break;
                }
            }

            if (IfEdgeChecked)
                continue;
        };

        int o1 = cuDFNsys::OrientationThree2DPnts(p1, q1, p2, _TOL_);
        int o2 = cuDFNsys::OrientationThree2DPnts(p1, q1, q2, _TOL_);
        int o3 = cuDFNsys::OrientationThree2DPnts(p2, q2, p1, _TOL_);
        int o4 = cuDFNsys::OrientationThree2DPnts(p2, q2, q1, _TOL_);

        if (stepNO == 1)
        {
            printf("particleNO: %d, edgeNO: %d, o: %d %d %d %d\n", particleNO, i, o1, o2, o3, o4);
        }

        if (o1 == 0 && o2 == 0 && o3 == 0 && o4 == 0)
        {
            // overlapping
            continue;
        };

        //float2 Edge2_[2] = {p2, q2};

        ///-------------------first end of the trajectory on an edge
        ///-------------------first end of the trajectory on an edge
        ///-------------------first end of the trajectory on an edge
        if (o3 == 0 && cuDFNsys::If2DPntLiesOnCollinearSeg(p2, p1, q2))
        {
            If_End1_on_Bound = true;
            edge_onBound_end1 = i;
            continue;
        }

        ///-------------------second end of the trajectory on an edge
        ///-------------------second end of the trajectory on an edge
        ///-------------------second end of the trajectory on an edge
        if (o4 == 0 && cuDFNsys::If2DPntLiesOnCollinearSeg(p2, q1, q2))
        {

            Result.z = (float)i;
            Result.x = q1.x;
            Result.y = q1.y;

            return Result;
        }

        // general case
        if (o1 != o2 && o3 != o4)
        {
            float norm_trajectory = sqrt((P_Trajectory[0].x - P_Trajectory[1].x) * (P_Trajectory[0].x - P_Trajectory[1].x) +
                                         (P_Trajectory[0].y - P_Trajectory[1].y) * (P_Trajectory[0].y - P_Trajectory[1].y));
            float norm_edge = sqrt((Vertex_Triangle[(i + 1) % 3].y - Vertex_Triangle[i].y) * (Vertex_Triangle[(i + 1) % 3].y - Vertex_Triangle[i].y) +
                                   (Vertex_Triangle[i].x - Vertex_Triangle[(i + 1) % 3].x) * (Vertex_Triangle[i].x - Vertex_Triangle[(i + 1) % 3].x));
            float factor_ = norm_edge / norm_trajectory;

            cuDFNsys::Scale2DSegment(P_Trajectory, factor_);

            float a1 = P_Trajectory[1].y - P_Trajectory[0].y;
            float b1 = P_Trajectory[0].x - P_Trajectory[1].x;
            float c1 = a1 * (P_Trajectory[0].x) + b1 * (P_Trajectory[0].y);

            // Line CD represented as a2x + b2y = c2
            float a2 = Vertex_Triangle[(i + 1) % 3].y - Vertex_Triangle[i].y;
            float b2 = Vertex_Triangle[i].x - Vertex_Triangle[(i + 1) % 3].x;
            float c2 = a2 * (Vertex_Triangle[i].x) + b2 * (Vertex_Triangle[i].y);

            float determinant = a1 * b2 - a2 * b1;

            Result.x = (b2 * c1 - b1 * c2) / determinant;
            Result.y = (a1 * c2 - a2 * c1) / determinant;

            Result.z = (float)i;

            return Result;
        };
    };

    // if end 1 is on bound, end 2 is out of bound
    if (If_End1_on_Bound == true)
    {
        Result.z = (float)edge_onBound_end1;
        Result.x = p1.x;
        Result.y = p1.y;

        return Result;
    }

    Result.z = -1.0f;
    return Result;
};
}; // namespace cuDFNsys