/****************************************************************************
* cuDFNsys - simulating flow and transport in 3D fracture networks          *
* Copyright (C) 2022, Tingchang YIN, Sergio GALINDO-TORRES                  *
*                                                                           *
* This program is free software: you can redistribute it and/or modify      *
* it under the terms of the GNU Affero General Public License as            *
* published by the Free Software Foundation, either version 3 of the        *
* License, or (at your option) any later version.                           *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU Affero General Public License for more details.                       *
*                                                                           *
* You should have received a copy of the GNU Affero General Public License  *
* along with this program.  If not, see <https://www.gnu.org/licenses/>.    *
*****************************************************************************/

#include "ParticleTransport/IdentifyParticleCrossesWhichEdge.cuh"

// ====================================================
// NAME:        IdentifyParticleCrossesWhichEdge
// DESCRIPTION: Identify which edge does the particle
//              cross
// AUTHOR:      Tingchang YIN
// DATE:        15/11/2022
// ====================================================

template <typename T>
__host__ __device__ cuDFNsys::Vector3<T> cuDFNsys::IdentifyParticleCrossesWhichEdge(cuDFNsys::Vector2<T> *P_Trajectory,
                                                                                    cuDFNsys::Vector2<T> *Vertex_Triangle,
                                                                                    T _TOL_,
                                                                                    cuDFNsys::Vector2<T> CrossedGlobalEdge[10][2],
                                                                                    int CountCrossedGlobalEdge,
                                                                                    uint stepNO,
                                                                                    uint particleNO)
{
    bool If_End1_on_Bound = false;
    int edge_onBound_end1 = -1;

    cuDFNsys::Vector3<T> Result;

    cuDFNsys::Vector2<T> p1, q1;
    p1 = P_Trajectory[0];
    q1 = P_Trajectory[1];

    for (uint i = 0; i < 3; ++i)
    {
        cuDFNsys::Vector2<T> p2, q2;
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
                cuDFNsys::Vector2<T> a = CrossedGlobalEdge[k][0];
                cuDFNsys::Vector2<T> b = CrossedGlobalEdge[k][1];

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

        int o1 = cuDFNsys::OrientationThree2DPnts<T>(p1, q1, p2, _TOL_);
        int o2 = cuDFNsys::OrientationThree2DPnts<T>(p1, q1, q2, _TOL_);
        int o3 = cuDFNsys::OrientationThree2DPnts<T>(p2, q2, p1, _TOL_);
        int o4 = cuDFNsys::OrientationThree2DPnts<T>(p2, q2, q1, _TOL_);

        //if (stepNO == 2000)
        //zprintf("edgeIDlocal: %d, o: %d %d %d %d\n", i, o1, o2, o3, o4);

        if (o1 == 0 && o2 == 0 && o3 == 0 && o4 == 0)
        {
            // overlapping
            continue;
        };

        //cuDFNsys::Vector2<T> Edge2_[2] = {p2, q2};

        ///-------------------first end of the trajectory on an edge
        ///-------------------first end of the trajectory on an edge
        ///-------------------first end of the trajectory on an edge
        if (o3 == 0 && cuDFNsys::If2DPntLiesOnCollinearSeg<T>(p2, p1, q2))
        {
            If_End1_on_Bound = true;
            edge_onBound_end1 = i;
            continue;
        }

        ///-------------------second end of the trajectory on an edge
        ///-------------------second end of the trajectory on an edge
        ///-------------------second end of the trajectory on an edge
        if (o4 == 0 && cuDFNsys::If2DPntLiesOnCollinearSeg<T>(p2, q1, q2))
        {

            Result.z = (T)i;
            Result.x = q1.x;
            Result.y = q1.y;

            return Result;
        }

        // general case
        if (o1 != o2 && o3 != o4)
        {
            // printf("xxxxxxxxxxx\n");
            // printf("LineSeg1:\n%.40f, %.40f\n%.40f, %.40f\n", P_Trajectory[0].x, P_Trajectory[0].y,
            //        P_Trajectory[1].x, P_Trajectory[1].y);
            // printf("LineSeg2:\n%.40f, %.40f\n%.40f, %.40f\n", Vertex_Triangle[i].x, Vertex_Triangle[i].y,
            //        Vertex_Triangle[(i + 1) % 3].x, Vertex_Triangle[(i + 1) % 3].y);

            // T norm_trajectory = sqrt((P_Trajectory[0].x - P_Trajectory[1].x) * (P_Trajectory[0].x - P_Trajectory[1].x) +
            //                          (P_Trajectory[0].y - P_Trajectory[1].y) * (P_Trajectory[0].y - P_Trajectory[1].y));
            // T norm_edge = sqrt((Vertex_Triangle[(i + 1) % 3].y - Vertex_Triangle[i].y) * (Vertex_Triangle[(i + 1) % 3].y - Vertex_Triangle[i].y) +
            //                    (Vertex_Triangle[i].x - Vertex_Triangle[(i + 1) % 3].x) * (Vertex_Triangle[i].x - Vertex_Triangle[(i + 1) % 3].x));
            //T factor_ = norm_edge / norm_trajectory;

            //cuDFNsys::Scale2DSegment<T>(P_Trajectory, factor_);

            // T a1 = P_Trajectory[1].y - P_Trajectory[0].y;
            // T b1 = P_Trajectory[0].x - P_Trajectory[1].x;
            // T c1 = a1 * (P_Trajectory[0].x) + b1 * (P_Trajectory[0].y);
            // // Line CD represented as a2x + b2y = c2
            // T a2 = Vertex_Triangle[(i + 1) % 3].y - Vertex_Triangle[i].y;
            // T b2 = Vertex_Triangle[i].x - Vertex_Triangle[(i + 1) % 3].x;
            // T c2 = a2 * (Vertex_Triangle[i].x) + b2 * (Vertex_Triangle[i].y);
            // T determinant = a1 * b2 - a2 * b1;
            // Result.x = (b2 * c1 - b1 * c2) / determinant;
            // Result.y = (a1 * c2 - a2 * c1) / determinant;
            // Result.z = (T)i;
            // printf("IdentifyParticleCrossesWhichEdge---Intersection: %.40f, %.40f\n", Result.x, Result.y);

            T F1 = P_Trajectory[1].x - Vertex_Triangle[(i + 1) % 3].x;
            T F2 = P_Trajectory[1].y - Vertex_Triangle[(i + 1) % 3].y;
            T M11 = P_Trajectory[1].x - P_Trajectory[0].x;
            T M21 = P_Trajectory[1].y - P_Trajectory[0].y;
            T M12 = Vertex_Triangle[i].x - Vertex_Triangle[(i + 1) % 3].x;
            T M22 = Vertex_Triangle[i].y - Vertex_Triangle[(i + 1) % 3].y;
            T deter = M11 * M22 - M12 * M21;
            T lambda = -(F2 * M12 - F1 * M22) / deter;
            Result.x = lambda * P_Trajectory[0].x + (1.0 - lambda) * P_Trajectory[1].x;
            Result.y = lambda * P_Trajectory[0].y + (1.0 - lambda) * P_Trajectory[1].y;
            Result.z = (T)i;
            // printf("IdentifyParticleCrossesWhichEdge---Intersection: %.40f, %.40f\n", Result.x, Result.y);
            return Result;
        };
    };

    //if (stepNO == 2000)
    //printf("-----------\n\n");

    // if end 1 is on bound, end 2 is out of bound
    if (If_End1_on_Bound == true)
    {
        Result.z = (T)edge_onBound_end1;
        Result.x = p1.x;
        Result.y = p1.y;

        return Result;
    }

    Result.z = -1.0f;
    return Result;
}; // IdentifyParticleCrossesWhichEdge
template __host__ __device__ cuDFNsys::Vector3<double> cuDFNsys::IdentifyParticleCrossesWhichEdge<double>(cuDFNsys::Vector2<double> *P_Trajectory,
                                                                                                          cuDFNsys::Vector2<double> *Vertex_Triangle,
                                                                                                          double _TOL_,
                                                                                                          cuDFNsys::Vector2<double> CrossedGlobalEdge[10][2],
                                                                                                          int CountCrossedGlobalEdge,
                                                                                                          uint stepNO,
                                                                                                          uint particleNO);
template __host__ __device__ cuDFNsys::Vector3<float> cuDFNsys::IdentifyParticleCrossesWhichEdge<float>(cuDFNsys::Vector2<float> *P_Trajectory,
                                                                                                        cuDFNsys::Vector2<float> *Vertex_Triangle,
                                                                                                        float _TOL_,
                                                                                                        cuDFNsys::Vector2<float> CrossedGlobalEdge[10][2],
                                                                                                        int CountCrossedGlobalEdge,
                                                                                                        uint stepNO,
                                                                                                        uint particleNO);