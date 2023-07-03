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

#include "ParticleTransport/ParticleMovementOneTimeStepGPUKernel.cuh"

// ====================================================
// NAME:        ParticleMovementOneTimeStepGPUKernel
// DESCRIPTION: ParticleMovementOneTimeStepGPUKernel
// AUTHOR:      Tingchang YIN
// DATE:        15/11/2022
// ====================================================
template <typename T>
__global__ void cuDFNsys::ParticleMovementOneTimeStepGPUKernel(unsigned long seed,
                                                               T delta_T_,
                                                               T Dispersion_local,
                                                               cuDFNsys::Particle<T> *P_DEV,
                                                               cuDFNsys::Fracture<T> *Frac_DEV,
                                                               cuDFNsys::EdgeToEle *EdgesSharedEle_DEV,
                                                               cuDFNsys::EleCoor<T> *Coordinate2D_Vec_dev_ptr,
                                                               //cuDFNsys::NeighborEle *NeighborEleOfOneEle_dev_ptr,
                                                               uint *EleToFracID_ptr,
                                                               T *velocity_ptr,
                                                               uint Dir_flow,
                                                               T outletcoordinate,
                                                               int count,
                                                               int numElements,
                                                               uint stepNO,
                                                               uint *Particle_runtime_error_dev_pnt,
                                                               uint NUMParticlesInTotal,
                                                               bool If_completeMixing)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= count)
        return;

    if (P_DEV[i].ParticleID < 0)
        return;

    uint EleID = P_DEV[i].ElementID; // from 1
    uint FracID;                     // from 0
    uint3 EdgeNO;                    // from 0
    cuDFNsys::Vector2<T> InitPos, TargPos;
    InitPos = P_DEV[i].Position2D;
    cuDFNsys::Vector3<T> Veloc_triangle;
    cuDFNsys::Vector2<T> Vertex_Triangle_ForVelocity[3];
    cuDFNsys::Vector2<T> Veloc_p;
    T conductivity_k;
    T b_aperture;

    //---------for calculated tortuosity
    cuDFNsys::Vector2<T> Pre_initPos = P_DEV[i].Position2D;
    uint Pre_fractureIDID = 0;
    //-----------------------------

    curandState state;

    curand_init(seed, i, 0, &state);

    T z1 = curand_normal(&state),
      z2 = curand_normal(&state);

    // if (i != 0)
    //     return;
    // // z1 = -3.7111978819675952578904798428993672132492;
    // // z2 = 3.4964017884433467031612963182851672172546;
    // InitPos.x = 0.6448155532070993789517387995147146284580;
    // InitPos.y = -1.5630180558642439159200421272544190287590;
    // EleID = 52528;

    bool IfUpdateTrajectoryLastStep = true;
    bool HaveRandomWalkerBeenThroughTrace = false;
    bool HaveRandomWalkerBeenReflected = false;

    T delta_T_consume = 0;
    T ratio_diffusion = 1;
    T diffusion_term_x = (T)1.0f * z1 * sqrt((T)2.0f * Dispersion_local * delta_T_),
      diffusion_term_y = (T)1.0f * z2 * sqrt((T)2.0f * Dispersion_local * delta_T_);

    T normLastTrajectory = 0;

    T _ParTran_MaxLoopTimes = 1000;

    for (uint Loop_time = 1;; Loop_time++)
    {
        if (Loop_time >= _ParTran_MaxLoopTimes)
        {
            printf("Particle %d, Loop times is too large!\n", i + 1);
            goto Debug100;
        }

        FracID = EleToFracID_ptr[EleID - 1];

        if (Loop_time == 1)
            Pre_fractureIDID = FracID;

        EdgeNO = make_uint3(EleID * 3 - 3, EleID * 3 - 2, EleID * 3 - 1);
        Veloc_triangle = cuDFNsys::MakeVector3(velocity_ptr[EdgeNO.x],
                                               velocity_ptr[EdgeNO.y],
                                               velocity_ptr[EdgeNO.z]);

        Vertex_Triangle_ForVelocity[0] = cuDFNsys::MakeVector2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[0], Coordinate2D_Vec_dev_ptr[EleID - 1].y[0]);
        Vertex_Triangle_ForVelocity[1] = cuDFNsys::MakeVector2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[1], Coordinate2D_Vec_dev_ptr[EleID - 1].y[1]);
        Vertex_Triangle_ForVelocity[2] = cuDFNsys::MakeVector2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[2], Coordinate2D_Vec_dev_ptr[EleID - 1].y[2]);

        if (IfUpdateTrajectoryLastStep)
        {
            Veloc_p = cuDFNsys::ReconstructVelocityGrid<T>(InitPos,
                                                           Vertex_Triangle_ForVelocity,
                                                           Veloc_triangle);
            conductivity_k = Frac_DEV[FracID].Conductivity;
            b_aperture = pow(conductivity_k * 12.0, 1.0 / 3.0);
            Veloc_p.x /= b_aperture, Veloc_p.y /= b_aperture;

            TargPos.x = InitPos.x + Veloc_p.x * (delta_T_ - delta_T_consume) + ratio_diffusion * diffusion_term_x;
            TargPos.y = InitPos.y + Veloc_p.y * (delta_T_ - delta_T_consume) + ratio_diffusion * diffusion_term_y;

            normLastTrajectory = pow((TargPos.x - InitPos.x) * (TargPos.x - InitPos.x) + (TargPos.y - InitPos.y) * (TargPos.y - InitPos.y), 0.5);
        }
        else
        {
            //IfUpdateTrajectoryLastStep = true;
        }

        //-------------debug---------
        {
            // printf("\nLoop_time: %d, Trajectory of particle ID %d:\n\t[%.40f, %.40f]\n\t[%.40f, %.40f]\n", Loop_time, i + 1,
            //        InitPos.x, InitPos.y, TargPos.x, TargPos.y);
            // printf("Loop_time: %d, Element ID %d of particle ID %d:\n\t[%.40f, %.40f]\n\t[%.40f, %.40f]\n\t[%.40f, %.40f]\n", Loop_time, EleID, i + 1,
            //        Vertex_Triangle_ForVelocity[0].x, Vertex_Triangle_ForVelocity[0].y, Vertex_Triangle_ForVelocity[1].x, Vertex_Triangle_ForVelocity[1].y,
            //        Vertex_Triangle_ForVelocity[2].x, Vertex_Triangle_ForVelocity[2].y);
            // cuDFNsys::Vector3<T> InitiPos3D = cuDFNsys::MakeVector3(InitPos.x, InitPos.y, (T)0.0);
            // cuDFNsys::Vector3<T> TargPos3D = cuDFNsys::MakeVector3(TargPos.x, TargPos.y, (T)0.0);
            // T RK_1[3][3];
            // Frac_DEV[FracID].RoationMatrix(RK_1, 23);
            // InitiPos3D = cuDFNsys::ProductSquare3Vector3<T>(RK_1, InitiPos3D);
            // TargPos3D = cuDFNsys::ProductSquare3Vector3<T>(RK_1, TargPos3D);
            // InitiPos3D.x += Frac_DEV[FracID].Center.x;
            // InitiPos3D.y += Frac_DEV[FracID].Center.y;
            // InitiPos3D.z += Frac_DEV[FracID].Center.z;
            // TargPos3D.x += Frac_DEV[FracID].Center.x;
            // TargPos3D.y += Frac_DEV[FracID].Center.y;
            // TargPos3D.z += Frac_DEV[FracID].Center.z;
            // printf("Trajectory3D of ParticleID %d, FracID: %d:\n\t%.40f, %.40f, %.40f\n\t%.40f, %.40f, %.40f\n",
            //        i + 1, FracID,
            //        InitiPos3D.x,
            //        InitiPos3D.y,
            //        InitiPos3D.z,
            //        TargPos3D.x,
            //        TargPos3D.y,
            //        TargPos3D.z);
        }
        //---------------------------

        bool IfTargPosInGrid = cuDFNsys::IfPntInside2DConvexPoly<T>(TargPos, Vertex_Triangle_ForVelocity, 3);
        // printf("IfTargPosInGrid: %d\n", IfTargPosInGrid);
        if (IfTargPosInGrid)
        {
            P_DEV[i].ElementID = EleID;
            P_DEV[i].Position2D = TargPos;
            // P_DEV[i].AccumDisplacement += pow((TargPos.x - InitPos.x) * (TargPos.x - InitPos.x) + (TargPos.y - InitPos.y) * (TargPos.y - InitPos.y), 0.5);
            goto TimeStepInfo;
        }

        uint2 IfOnBorder = cuDFNsys::IfRandomWalkLieOnBorderTriangle<T>(TargPos, Vertex_Triangle_ForVelocity, (T)0.0);

        if (IfOnBorder.x == 1)
        {
            // printf("Error: the trajectory (particle ID: %d) starts at the vertex\n", i + 1);
            // printf("Trajectory of particle ID %d:\n\t[%.40f, %.40f]\n\t[%.40f, %.40f]\n", i + 1,
            //        InitPos.x, InitPos.y, TargPos.x, TargPos.y);
            // printf("Element ID %d of particle ID %d:\n\t[%.40f, %.40f]\n\t[%.40f, %.40f]\n\t[%.40f, %.40f]\n", EleID, i + 1,
            //        Vertex_Triangle_ForVelocity[0].x, Vertex_Triangle_ForVelocity[0].y, Vertex_Triangle_ForVelocity[1].x, Vertex_Triangle_ForVelocity[1].y,
            //        Vertex_Triangle_ForVelocity[2].x, Vertex_Triangle_ForVelocity[2].y);
            // goto Debug100;
            // check if reached

            P_DEV[i].ElementID = EleID;
            P_DEV[i].Position2D = TargPos;

            // if this target point reaches the outlet plane?
            cuDFNsys::Vector3<T> Pos3D_ = cuDFNsys::Roate2DPositionTo3D<T>(TargPos,
                                                                           Frac_DEV[EleToFracID_ptr[EleID - 1]]);
            T *tmp_pnt = &(Pos3D_.x);
            if (abs(tmp_pnt[Dir_flow] - outletcoordinate) < 1e-4)
            {
                P_DEV[i].ParticleID = -P_DEV[i].ParticleID; // reaches outlet plane
                // P_DEV[i].AccumDisplacement += pow((TargPos.x - InitPos.x) * (TargPos.x - InitPos.x) + (TargPos.y - InitPos.y) * (TargPos.y - InitPos.y), 0.5);

                goto TimeStepInfo;
            }
            else
            {
                // P_DEV[i].AccumDisplacement += pow((TargPos.x - InitPos.x) * (TargPos.x - InitPos.x) + (TargPos.y - InitPos.y) * (TargPos.y - InitPos.y), 0.5);
                // just on bounds of grids
                goto TimeStepInfo;
            }
        }

        //---------------outside
        bool IFintersectBorder;
        cuDFNsys::Vector2<T> Trajectory[2] = {InitPos, TargPos};

        // rescale the trajectory
        Trajectory[0].x = InitPos.x + (TargPos.x - InitPos.x) * 0.001;
        Trajectory[0].y = InitPos.y + (TargPos.y - InitPos.y) * 0.001;

        uint result4[4] = {0};
        IFintersectBorder = cuDFNsys::WhichEdgeDoesTrajectoryIntersect<T>(Trajectory,
                                                                          Vertex_Triangle_ForVelocity,
                                                                          result4);

        // printf("rescaled Trajectory of particle ID %d:\n\t[%.40f, %.40f]\n\t[%.40f, %.40f]\n", i + 1,
        // Trajectory[0].x, Trajectory[0].y, TargPos.x, TargPos.y);
        if (!IFintersectBorder)
        {
            uint MaxLoopInterect = 30;
            for (uint yy = 1; yy <= MaxLoopInterect; ++yy)
            {
                T U_ratio = (yy == MaxLoopInterect ? 0 : 0.001 / (pow(10.0, (T)yy)));

                Trajectory[0].x = InitPos.x + (TargPos.x - InitPos.x) * U_ratio;
                Trajectory[0].y = InitPos.y + (TargPos.y - InitPos.y) * U_ratio;
                //Trajectory[1] = TargPos;
                IFintersectBorder = cuDFNsys::WhichEdgeDoesTrajectoryIntersect<T>(Trajectory,
                                                                                  Vertex_Triangle_ForVelocity,
                                                                                  result4);

                //printf("y: %d, rescaled Trajectory of particle ID %d:\n\t[%.40f, %.40f]\n\t[%.40f, %.40f]\n", yy, i + 1,
                // Trajectory[0].x, Trajectory[0].y, TargPos.x, TargPos.y);

                if (IFintersectBorder && result4[0] == 1)
                    break;
            }

            if (!IFintersectBorder)
            {
                // the initial point may lies on the border of the triangle
                // and the target position is outside the triangle
                // but the algorithm cannot determine the which edge does the trajectory intersect with
                // so , elongate the trajectory a little bit
                Trajectory[0].x = InitPos.x - (TargPos.x - InitPos.x) * 0.0000001;
                Trajectory[0].y = InitPos.y - (TargPos.y - InitPos.y) * 0.0000001;
                IFintersectBorder = cuDFNsys::WhichEdgeDoesTrajectoryIntersect<T>(Trajectory,
                                                                                  Vertex_Triangle_ForVelocity,
                                                                                  result4);
                if (!IFintersectBorder)
                {
                    printf("Error: did not find which edge of the triangle that intersects with the trajectory (particle ID: %d)!\n", i + 1);
                    printf("Trajectory of particle ID %d:\n\t[%.40f, %.40f]\n\t[%.40f, %.40f]\n", i + 1,
                           InitPos.x, InitPos.y, TargPos.x, TargPos.y);
                    printf("Element ID %d of particle ID %d:\n\t[%.40f, %.40f]\n\t[%.40f, %.40f]\n\t[%.40f, %.40f]\n", EleID, i + 1,
                           Vertex_Triangle_ForVelocity[0].x, Vertex_Triangle_ForVelocity[0].y, Vertex_Triangle_ForVelocity[1].x, Vertex_Triangle_ForVelocity[1].y,
                           Vertex_Triangle_ForVelocity[2].x, Vertex_Triangle_ForVelocity[2].y);
                    goto Debug100;
                }
            }
        }

        if (result4[0] > 1)
        {
            if (IFintersectBorder && result4[0] == 2)
            {
                // the trajectory exactly across a vertex
                // let the initial position be the vertex
                // address it at border
                // T norm_FinishedTrajec = pow((Vertex_Triangle_ForVelocity[result4[2]].x - InitPos.x) * (Vertex_Triangle_ForVelocity[result4[2]].x - InitPos.x) +
                //                                 (Vertex_Triangle_ForVelocity[result4[2]].y - InitPos.y) *
                //                                     (Vertex_Triangle_ForVelocity[result4[2]].y - InitPos.y),
                //                             0.5);
                // T ratio_Trajectory = norm_FinishedTrajec / normLastTrajectory;
                // delta_T_consume += (delta_T_ - delta_T_consume) * ratio_Trajectory;
                // ratio_diffusion = (1.0 - delta_T_consume / delta_T_);
                // InitPos = Vertex_Triangle_ForVelocity[result4[2]];
                // P_DEV[i].AccumDisplacement += norm_FinishedTrajec;
                // IfStartFromVertexAndOutsideElement = true;
                // continue;
                printf("Error: the trajectory (particle ID: %d) may accross a vertex of element %d? Impossible!\n", i + 1, EleID);
                printf("Trajectory of particle ID %d:\n\t[%.40f, %.40f]\n\t[%.40f, %.40f]\n", i + 1,
                       InitPos.x, InitPos.y, TargPos.x, TargPos.y);
                printf("Element ID %d of particle ID %d:\n\t[%.40f, %.40f]\n\t[%.40f, %.40f]\n\t[%.40f, %.40f]\n", EleID, i + 1,
                       Vertex_Triangle_ForVelocity[0].x, Vertex_Triangle_ForVelocity[0].y, Vertex_Triangle_ForVelocity[1].x, Vertex_Triangle_ForVelocity[1].y,
                       Vertex_Triangle_ForVelocity[2].x, Vertex_Triangle_ForVelocity[2].y);
                goto Debug100;
            }
            else if (IFintersectBorder && result4[0] > 2)
            {
                printf("Error: the trajectory (particle ID: %d) intersects to three edges of a triangle? Impossible!\n", i + 1);
                printf("Trajectory of particle ID %d:\n\t[%.40f, %.40f]\n\t[%.40f, %.40f]\n", i + 1,
                       InitPos.x, InitPos.y, TargPos.x, TargPos.y);
                printf("Element ID %d of particle ID %d:\n\t[%.40f, %.40f]\n\t[%.40f, %.40f]\n\t[%.40f, %.40f]\n", EleID, i + 1,
                       Vertex_Triangle_ForVelocity[0].x, Vertex_Triangle_ForVelocity[0].y, Vertex_Triangle_ForVelocity[1].x, Vertex_Triangle_ForVelocity[1].y,
                       Vertex_Triangle_ForVelocity[2].x, Vertex_Triangle_ForVelocity[2].y);
                goto Debug100;
            }
        }

        cuDFNsys::Vector2<T> TheEdge[2] = {Vertex_Triangle_ForVelocity[result4[1]],
                                           Vertex_Triangle_ForVelocity[(result4[1] + 1) % 3]};

        cuDFNsys::Vector3<T> Inj = cuDFNsys::IntersectionLineLine2D<T>(TheEdge, Trajectory);

        if (Inj.x == 0)
        {
            printf("Error: the trajectory (particle ID: %d) is not intersecting with the edge?\n", i + 1);
            printf("Trajectory of particle ID %d:\n\t[%.40f, %.40f]\n\t[%.40f, %.40f]\n", i + 1,
                   InitPos.x, InitPos.y, TargPos.x, TargPos.y);
            printf("The edge of particle ID %d intersected:\n\t[%.40f, %.40f]\n\t[%.40f, %.40f]\n", i + 1,
                   TheEdge[0].x, TheEdge[0].y, TheEdge[1].x, TheEdge[1].y);
            printf("Element ID %d of particle ID %d:\n\t[%.40f, %.40f]\n\t[%.40f, %.40f]\n\t[%.40f, %.40f]\n", EleID, i + 1,
                   Vertex_Triangle_ForVelocity[0].x, Vertex_Triangle_ForVelocity[0].y, Vertex_Triangle_ForVelocity[1].x, Vertex_Triangle_ForVelocity[1].y,
                   Vertex_Triangle_ForVelocity[2].x, Vertex_Triangle_ForVelocity[2].y);
            goto Debug100;
        }

        cuDFNsys::Vector2<T> IntersectionOnEdge = cuDFNsys::MakeVector2(Inj.y, Inj.z);

        uint *tyu = &(EdgeNO.x);
        uint GlobalEdgeNO = tyu[result4[1]];

        uint NumOfElesSharedEdge = EdgesSharedEle_DEV[GlobalEdgeNO].NumSharedEle;

        if (NumOfElesSharedEdge == 1)
        {
            cuDFNsys::Vector3<T> Pos3D_ = cuDFNsys::Roate2DPositionTo3D<T>(IntersectionOnEdge,
                                                                           Frac_DEV[EleToFracID_ptr[EleID - 1]]);
            T *tmp_pnt = &(Pos3D_.x);

            if (abs(tmp_pnt[Dir_flow] - outletcoordinate) < 1e-4)
            {
                P_DEV[i].ParticleID = -P_DEV[i].ParticleID; // reaches outlet plane
                // P_DEV[i].AccumDisplacement += pow((IntersectionOnEdge.x - InitPos.x) * (IntersectionOnEdge.x - InitPos.x) + (IntersectionOnEdge.y - InitPos.y) * (IntersectionOnEdge.y - InitPos.y), 0.5);

                goto TimeStepInfo;
            }
            else if (tmp_pnt[Dir_flow] > ((T)-1.0 * outletcoordinate - 1e-4))
            {
                // go out from the inlet
                if (Dispersion_local != 0)
                {
                    P_DEV[i].ParticleID = NUMParticlesInTotal * 2; //delete this randam walker
                    goto TimeStepInfo;
                }
                else
                {
                    if (!HaveRandomWalkerBeenThroughTrace)
                    {
                        printf("Error: the trajectory (particle ID: %d) is acrossing the inlet, but the molecular diffusion is zero\n", i + 1);
                        printf("Trajectory of particle ID %d:\n\t[%.40f, %.40f]\n\t[%.40f, %.40f]\n", i + 1,
                               InitPos.x, InitPos.y, TargPos.x, TargPos.y);
                        printf("The edge of particle ID %d intersected:\n\t[%.40f, %.40f]\n\t[%.40f, %.40f]\n", i + 1,
                               TheEdge[0].x, TheEdge[0].y, TheEdge[1].x, TheEdge[1].y);
                        printf("Intersection of particle ID %d intersected:\n\t[%.40f, %.40f]\n", i + 1,
                               IntersectionOnEdge.x, IntersectionOnEdge.y);
                        printf("Element ID %d of particle ID %d:\n\t[%.40f, %.40f]\n\t[%.40f, %.40f]\n\t[%.40f, %.40f]\n", EleID, i + 1,
                               Vertex_Triangle_ForVelocity[0].x, Vertex_Triangle_ForVelocity[0].y, Vertex_Triangle_ForVelocity[1].x, Vertex_Triangle_ForVelocity[1].y,
                               Vertex_Triangle_ForVelocity[2].x, Vertex_Triangle_ForVelocity[2].y);
                        goto Debug100;
                    }
                    else
                    {
                        T norm_FinishedTrajec = pow((IntersectionOnEdge.x - InitPos.x) * (IntersectionOnEdge.x - InitPos.x) + (IntersectionOnEdge.y - InitPos.y) * (IntersectionOnEdge.y - InitPos.y), 0.5);
                        T ratio_Trajectory = norm_FinishedTrajec / normLastTrajectory;
                        delta_T_consume += (delta_T_ - delta_T_consume) * ratio_Trajectory;
                        ratio_diffusion = (1.0 - delta_T_consume / delta_T_);

                        // P_DEV[i].AccumDisplacement += norm_FinishedTrajec;

                        InitPos = IntersectionOnEdge;
                        IfUpdateTrajectoryLastStep = true;

                        continue;
                    }
                }
            }
            else // non-flux edge
            {
                T norm_FinishedTrajec = pow((IntersectionOnEdge.x - InitPos.x) * (IntersectionOnEdge.x - InitPos.x) + (IntersectionOnEdge.y - InitPos.y) * (IntersectionOnEdge.y - InitPos.y), 0.5);
                T ratio_Trajectory = norm_FinishedTrajec / normLastTrajectory;
                delta_T_consume += (delta_T_ - delta_T_consume) * ratio_Trajectory;
                ratio_diffusion = (1.0 - delta_T_consume / delta_T_);

                // P_DEV[i].AccumDisplacement += norm_FinishedTrajec;

                // printf("Reflect at non-flux boundary, particle ID %d: \n", i + 1);
                // printf("The edge of particle ID %d intersected:\n\t[%.40f, %.40f]\n\t[%.40f, %.40f]\n", i + 1,
                //        TheEdge[0].x, TheEdge[0].y, TheEdge[1].x, TheEdge[1].y);
                // printf("Intersection of particle ID %d intersected:\n\t[%.40f, %.40f]\n", i + 1,
                //        IntersectionOnEdge.x, IntersectionOnEdge.y);
                // printf("Element ID %d of particle ID %d:\n\t[%.40f, %.40f]\n\t[%.40f, %.40f]\n\t[%.40f, %.40f]\n", EleID, i + 1,
                //        Vertex_Triangle_ForVelocity[0].x, Vertex_Triangle_ForVelocity[0].y, Vertex_Triangle_ForVelocity[1].x, Vertex_Triangle_ForVelocity[1].y,
                //        Vertex_Triangle_ForVelocity[2].x, Vertex_Triangle_ForVelocity[2].y);
                // printf("Trajectory of particle ID %d:\n\t[%.40f, %.40f]\n\t[%.40f, %.40f]\n", i + 1,
                //        InitPos.x, InitPos.y, TargPos.x, TargPos.y);

                //----------------
                TargPos = cuDFNsys::ParticleReflection<T>(TargPos, Vertex_Triangle_ForVelocity[result4[1]],
                                                          Vertex_Triangle_ForVelocity[(result4[1] + 1) % 3]);
                InitPos = IntersectionOnEdge;

                // printf("after the reflection, Trajectory of particle ID %d:\n\t[%.40f, %.40f]\n\t[%.40f, %.40f]\n", i + 1,
                //        InitPos.x, InitPos.y, TargPos.x, TargPos.y);

                IfUpdateTrajectoryLastStep = false;
                normLastTrajectory = pow((TargPos.x - InitPos.x) * (TargPos.x - InitPos.x) + (TargPos.y - InitPos.y) * (TargPos.y - InitPos.y), 0.5);
                // goto Debug100;
                HaveRandomWalkerBeenReflected = true;

                continue;
            }
        }
        else if (NumOfElesSharedEdge == 2)
        {
            T norm_FinishedTrajec = pow((IntersectionOnEdge.x - InitPos.x) * (IntersectionOnEdge.x - InitPos.x) + (IntersectionOnEdge.y - InitPos.y) * (IntersectionOnEdge.y - InitPos.y), 0.5);
            T ratio_Trajectory = norm_FinishedTrajec / normLastTrajectory;
            delta_T_consume += (delta_T_ - delta_T_consume) * ratio_Trajectory;
            ratio_diffusion = (1.0 - delta_T_consume / delta_T_);

            InitPos = IntersectionOnEdge;
            EleID = (EdgesSharedEle_DEV[GlobalEdgeNO].EleID[0] == EleID ? EdgesSharedEle_DEV[GlobalEdgeNO].EleID[1] : EdgesSharedEle_DEV[GlobalEdgeNO].EleID[0]);

            // P_DEV[i].AccumDisplacement += norm_FinishedTrajec;

            if ((HaveRandomWalkerBeenThroughTrace || HaveRandomWalkerBeenReflected) && Dispersion_local == 0)
                IfUpdateTrajectoryLastStep = true;

            continue;
        }
        else
        {
            // printf("NumOfElesSharedEdge of Particle ID %d is %d\n", i + 1, NumOfElesSharedEdge);
            //goto Debug100;
            // may go to another fracture
            // may go to another fracture
            // may go to another fracture

            uint PreviousFracID = FracID;
            //uint PreEleID = EleID;
            // determine where to go
            // determine where to go
            // determine where to go
            int IndexLocal = -1; // index of local shared edge NO the next element;
            // the edge is the one where the intersection lies on

            int newELEID_ = -1;
            int newFracID_ = -1;

            //for (uint p = 0; p < NumOfElesSharedEdge; ++p)
            //printf("EdgesSharedEle_DEV[GlobalEdgeNO].EleID[p]: %d, fracid: %d\n", EdgesSharedEle_DEV[GlobalEdgeNO].EleID[p], EleToFracID_ptr[EdgesSharedEle_DEV[GlobalEdgeNO].EleID[p] - 1]);

            bool ifAllsharedEdgeVelocityPositive = false;
            cuDFNsys::WhichElementToGo<T>(EleID,
                                          NumOfElesSharedEdge,
                                          Dispersion_local,
                                          EdgesSharedEle_DEV[GlobalEdgeNO].EleID,
                                          EdgesSharedEle_DEV[GlobalEdgeNO].LocalEdgeNO,
                                          EleToFracID_ptr,
                                          Frac_DEV,
                                          Coordinate2D_Vec_dev_ptr,
                                          velocity_ptr,
                                          curand_uniform(&state),
                                          newELEID_,
                                          newFracID_,
                                          IndexLocal,
                                          ifAllsharedEdgeVelocityPositive, If_completeMixing);
            // printf("ifAllsharedEdgeVelocityPositive: %d\n", ifAllsharedEdgeVelocityPositive);

            if (ifAllsharedEdgeVelocityPositive == false && (newELEID_ == -1 || newFracID_ == -1 || IndexLocal == -1))
            {
                printf("Warning: illeagal indice appear, because the determination of which element to go is failed!\nnewELEID_: %d\nnewFracID_: %d\nIndexLocal: %d\n",
                       newELEID_,
                       newFracID_,
                       IndexLocal);
                goto Debug100;
            }

            T norm_FinishedTrajec = pow((IntersectionOnEdge.x - InitPos.x) * (IntersectionOnEdge.x - InitPos.x) + (IntersectionOnEdge.y - InitPos.y) * (IntersectionOnEdge.y - InitPos.y), 0.5);
            T ratio_Trajectory = norm_FinishedTrajec / normLastTrajectory;
            delta_T_consume += (delta_T_ - delta_T_consume) * ratio_Trajectory;
            ratio_diffusion = (1.0 - delta_T_consume / delta_T_);

            if (ifAllsharedEdgeVelocityPositive == true)
            {
                /// the particle trajectory goes back to the previous element, by reflection
                /// the particle trajectory goes back to the previous element, by reflection
                /// the particle trajectory goes back to the previous element, by reflection
                //printf("\nifAllsharedEdgeVelocityPositive == true\n");
                InitPos = IntersectionOnEdge;
                TargPos = cuDFNsys::ParticleReflection<T>(TargPos, TheEdge[0], TheEdge[1]);

                IfUpdateTrajectoryLastStep = false;
                // P_DEV[i].AccumDisplacement += norm_FinishedTrajec;

                normLastTrajectory = pow((TargPos.x - InitPos.x) * (TargPos.x - InitPos.x) + (TargPos.y - InitPos.y) * (TargPos.y - InitPos.y), 0.5);
                HaveRandomWalkerBeenReflected = true;

                continue;
            };

            FracID = (uint)newFracID_;
            EleID = (uint)newELEID_;
            EdgeNO = make_uint3(EleID * 3 - 3, EleID * 3 - 2, EleID * 3 - 1); // from 0

            if (FracID == PreviousFracID) // did not change fracture ID
            {
                // InitPos = IntersectionOnEdge;

                InitPos = IntersectionOnEdge;
                // P_DEV[i].AccumDisplacement += norm_FinishedTrajec;

                continue;
            }

            //-------------new------------- intersection point
            T RK_1[3][3];
            Frac_DEV[PreviousFracID].RoationMatrix(RK_1, 23);
            cuDFNsys::Vector3<T> IntersectionOnEdge3D = cuDFNsys::MakeVector3(IntersectionOnEdge.x, IntersectionOnEdge.y, (T)0);
            IntersectionOnEdge3D = cuDFNsys::ProductSquare3Vector3<T>(RK_1, IntersectionOnEdge3D);
            IntersectionOnEdge3D = cuDFNsys::MakeVector3(IntersectionOnEdge3D.x + Frac_DEV[PreviousFracID].Center.x,
                                                         IntersectionOnEdge3D.y + Frac_DEV[PreviousFracID].Center.y,
                                                         IntersectionOnEdge3D.z + Frac_DEV[PreviousFracID].Center.z);

            cuDFNsys::Vector3<T> IntersectionOnEdge3D_anotherFrac = cuDFNsys::MakeVector3(IntersectionOnEdge3D.x - Frac_DEV[FracID].Center.x,
                                                                                          IntersectionOnEdge3D.y - Frac_DEV[FracID].Center.y,
                                                                                          IntersectionOnEdge3D.z - Frac_DEV[FracID].Center.z);

            T RK_2[3][3];
            Frac_DEV[FracID].RoationMatrix(RK_2, 32);
            IntersectionOnEdge3D_anotherFrac = cuDFNsys::ProductSquare3Vector3<T>(RK_2, IntersectionOnEdge3D_anotherFrac);
            cuDFNsys::Vector2<T> IntersectionOnEdge2D_anotherFrac = cuDFNsys::MakeVector2(IntersectionOnEdge3D_anotherFrac.x,
                                                                                          IntersectionOnEdge3D_anotherFrac.y);
            //---------------------------------------------------
            // TheEdge 3D
            cuDFNsys::Vector3<T> TheEdge_oldEle[2];
            TheEdge_oldEle[0] = cuDFNsys::MakeVector3(TheEdge[0].x, TheEdge[0].y, (T)0);
            TheEdge_oldEle[1] = cuDFNsys::MakeVector3(TheEdge[1].x, TheEdge[1].y, (T)0);
            TheEdge_oldEle[0] = cuDFNsys::ProductSquare3Vector3<T>(RK_1, TheEdge_oldEle[0]);
            TheEdge_oldEle[1] = cuDFNsys::ProductSquare3Vector3<T>(RK_1, TheEdge_oldEle[1]);
            TheEdge_oldEle[0] = cuDFNsys::MakeVector3(TheEdge_oldEle[0].x + Frac_DEV[PreviousFracID].Center.x,
                                                      TheEdge_oldEle[0].y + Frac_DEV[PreviousFracID].Center.y,
                                                      TheEdge_oldEle[0].z + Frac_DEV[PreviousFracID].Center.z);
            TheEdge_oldEle[1] = cuDFNsys::MakeVector3(TheEdge_oldEle[1].x + Frac_DEV[PreviousFracID].Center.x,
                                                      TheEdge_oldEle[1].y + Frac_DEV[PreviousFracID].Center.y,
                                                      TheEdge_oldEle[1].z + Frac_DEV[PreviousFracID].Center.z);
            //-------------------
            // new edge 3D
            cuDFNsys::Vector2<T> TriangleNew_ele[3];
            TriangleNew_ele[0] = cuDFNsys::MakeVector2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[0], Coordinate2D_Vec_dev_ptr[EleID - 1].y[0]);
            TriangleNew_ele[1] = cuDFNsys::MakeVector2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[1], Coordinate2D_Vec_dev_ptr[EleID - 1].y[1]);
            TriangleNew_ele[2] = cuDFNsys::MakeVector2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[2], Coordinate2D_Vec_dev_ptr[EleID - 1].y[2]);

            cuDFNsys::Vector3<T> TheEdge_newEle[2];
            TheEdge_newEle[0] = cuDFNsys::MakeVector3(TriangleNew_ele[IndexLocal].x, TriangleNew_ele[IndexLocal].y, (T)0);
            TheEdge_newEle[1] = cuDFNsys::MakeVector3(TriangleNew_ele[(IndexLocal + 1) % 3].x, TriangleNew_ele[(IndexLocal + 1) % 3].y, (T)0);

            T RK_3[3][3];
            Frac_DEV[FracID].RoationMatrix(RK_3, 23);
            TheEdge_newEle[0] = cuDFNsys::ProductSquare3Vector3<T>(RK_3, TheEdge_newEle[0]);
            TheEdge_newEle[1] = cuDFNsys::ProductSquare3Vector3<T>(RK_3, TheEdge_newEle[1]);
            TheEdge_newEle[0] = cuDFNsys::MakeVector3(TheEdge_newEle[0].x + Frac_DEV[FracID].Center.x,
                                                      TheEdge_newEle[0].y + Frac_DEV[FracID].Center.y,
                                                      TheEdge_newEle[0].z + Frac_DEV[FracID].Center.z);
            TheEdge_newEle[1] = cuDFNsys::MakeVector3(TheEdge_newEle[1].x + Frac_DEV[FracID].Center.x,
                                                      TheEdge_newEle[1].y + Frac_DEV[FracID].Center.y,
                                                      TheEdge_newEle[1].z + Frac_DEV[FracID].Center.z);

            //----------------
            // if two edges are oppsite?
            bool IfTwoEdgesOpposite = false;
            if (abs(sqrt((TheEdge_oldEle[0].x - TheEdge_newEle[0].x) * (TheEdge_oldEle[0].x - TheEdge_newEle[0].x) + (TheEdge_oldEle[0].y - TheEdge_newEle[0].y) * (TheEdge_oldEle[0].y - TheEdge_newEle[0].y))) < 1e-3 &&
                abs(sqrt((TheEdge_oldEle[1].x - TheEdge_newEle[1].x) * (TheEdge_oldEle[1].x - TheEdge_newEle[1].x) + (TheEdge_oldEle[1].y - TheEdge_newEle[1].y) * (TheEdge_oldEle[1].y - TheEdge_newEle[1].y))) < 1e-3)
            {
                IfTwoEdgesOpposite = false;
            }
            else if (abs(sqrt((TheEdge_oldEle[0].x - TheEdge_newEle[1].x) * (TheEdge_oldEle[0].x - TheEdge_newEle[1].x) + (TheEdge_oldEle[0].y - TheEdge_newEle[1].y) * (TheEdge_oldEle[0].y - TheEdge_newEle[1].y))) < 1e-3 &&
                     abs(sqrt((TheEdge_oldEle[1].x - TheEdge_newEle[0].x) * (TheEdge_oldEle[1].x - TheEdge_newEle[0].x) + (TheEdge_oldEle[1].y - TheEdge_newEle[0].y) * (TheEdge_oldEle[1].y - TheEdge_newEle[0].y))) < 1e-3)
            {
                IfTwoEdgesOpposite = true;
            }
            else
            {
                printf("Error: particle ID: %d, when a random walker get through a trace, the shared edge is not accurately identified\n", i + 1);
                goto Debug100;
            }

            // the new edge in 2D
            cuDFNsys::Vector2<T> TheEdge_new[2] = {cuDFNsys::MakeVector2(TriangleNew_ele[IndexLocal].x, TriangleNew_ele[IndexLocal].y),
                                                   cuDFNsys::MakeVector2(TriangleNew_ele[(IndexLocal + 1) % 3].x, TriangleNew_ele[(IndexLocal + 1) % 3].y)};

            if (!IfTwoEdgesOpposite)
            {
                TheEdge_new[1] = cuDFNsys::MakeVector2(TriangleNew_ele[IndexLocal].x, TriangleNew_ele[IndexLocal].y);
                TheEdge_new[0] = cuDFNsys::MakeVector2(TriangleNew_ele[(IndexLocal + 1) % 3].x, TriangleNew_ele[(IndexLocal + 1) % 3].y);
            }

            // the angle between the old edge and the old remainning trajectory
            cuDFNsys::Vector2<T> vector_1_edge = cuDFNsys::MakeVector2(TheEdge[1].x - TheEdge[0].x, TheEdge[1].y - TheEdge[0].y),
                                 vector_1_trajectOld = cuDFNsys::MakeVector2(TargPos.x - IntersectionOnEdge.x, TargPos.y - IntersectionOnEdge.y);
            T norm_1_edge = sqrt(vector_1_edge.x * vector_1_edge.x + vector_1_edge.y * vector_1_edge.y);
            T norm_2_trjOld = sqrt(vector_1_trajectOld.x * vector_1_trajectOld.x + vector_1_trajectOld.y * vector_1_trajectOld.y);
            vector_1_edge.x /= norm_1_edge, vector_1_edge.y /= norm_1_edge;
            vector_1_trajectOld.x /= norm_2_trjOld, vector_1_trajectOld.y /= norm_2_trjOld;
            // https://wumbo.net/formulas/angle-between-two-vectors-2d/
            // angle amount from vector 1 to vector 2
            T angle_s = atan2(vector_1_trajectOld.y * vector_1_edge.x - vector_1_trajectOld.x * vector_1_edge.y,
                              vector_1_trajectOld.x * vector_1_edge.x + vector_1_trajectOld.y * vector_1_edge.y);

            // map it to the new element
            cuDFNsys::Vector2<T> vector_2_edge = cuDFNsys::MakeVector2(TheEdge_new[1].x - TheEdge_new[0].x, TheEdge_new[1].y - TheEdge_new[0].y);
            T norma_2_edge = sqrt(vector_2_edge.x * vector_2_edge.x + vector_2_edge.y * vector_2_edge.y);
            vector_2_edge.x /= norma_2_edge, vector_2_edge.y /= norma_2_edge;

            bool IfHasBeenReverse = false;
        reverseDir:;

            cuDFNsys::Quaternion<T> qua;
            qua = qua.DescribeRotation(cuDFNsys::MakeVector3((T)0, (T)0, (T)1), angle_s);
            cuDFNsys::Vector3<T> temp_new_tra = qua.Rotate(cuDFNsys::MakeVector3(vector_2_edge.x, vector_2_edge.y, (T)0));

            cuDFNsys::Vector2<T> newTarget = cuDFNsys::MakeVector2(temp_new_tra.x, temp_new_tra.y);
            T norm_new_targ = sqrt(newTarget.x * newTarget.x + newTarget.y * newTarget.y);
            newTarget.x /= norm_new_targ, newTarget.y /= norm_new_targ;

            //newTarget.x += IntersectionOnEdge2D_anotherFrac.x;
            //newTarget.y += IntersectionOnEdge2D_anotherFrac.y;

            newTarget.x = IntersectionOnEdge2D_anotherFrac.x + newTarget.x * norm_2_trjOld;
            newTarget.y = IntersectionOnEdge2D_anotherFrac.y + newTarget.y * norm_2_trjOld;

            InitPos = IntersectionOnEdge2D_anotherFrac;
            TargPos = newTarget;

            // a check; a check; a check; a check; a check.
            if (!IfHasBeenReverse)
            {
                // vector_2_edge
                cuDFNsys::Vector2<T> vector_2_ase = cuDFNsys::MakeVector2(TargPos.x - InitPos.x,
                                                                          TargPos.y - InitPos.y);
                T norm_vector_2_ase = sqrt(vector_2_ase.x * vector_2_ase.x + vector_2_ase.y * vector_2_ase.y);
                vector_2_ase.x /= norm_vector_2_ase, vector_2_ase.y /= norm_vector_2_ase;

                // one
                cuDFNsys::Vector2<T> vector_3_ase = cuDFNsys::MakeVector2(TriangleNew_ele[(IndexLocal + 2) % 3].x - TheEdge_new[0].x,
                                                                          TriangleNew_ele[(IndexLocal + 2) % 3].y - TheEdge_new[0].y);
                norm_vector_2_ase = sqrt(vector_3_ase.x * vector_3_ase.x + vector_3_ase.y * vector_3_ase.y);
                vector_3_ase.x /= norm_vector_2_ase, vector_3_ase.y /= norm_vector_2_ase;

                // angle
                T angle_u1 = atan2(vector_2_ase.y * vector_2_edge.x - vector_2_ase.x * vector_2_edge.y,
                                   vector_2_ase.x * vector_2_edge.x + vector_2_ase.y * vector_2_edge.y);
                T angle_u2 = atan2(vector_3_ase.y * vector_2_edge.x - vector_3_ase.x * vector_2_edge.y,
                                   vector_3_ase.x * vector_2_edge.x + vector_3_ase.y * vector_2_edge.y);
                if (cuDFNsys::Sgn<T>(angle_u1) != cuDFNsys::Sgn<T>(angle_u2) && !IfHasBeenReverse)
                {
                    angle_s = -angle_s;
                    IfHasBeenReverse = true;
                    goto reverseDir;
                    /// Vertex_Triangle_ForVelocity[0] = cuDFNsys::MakeVector2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[0], Coordinate2D_Vec_dev_ptr[EleID - 1].y[0]);
                    /// Vertex_Triangle_ForVelocity[1] = cuDFNsys::MakeVector2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[1], Coordinate2D_Vec_dev_ptr[EleID - 1].y[1]);
                    /// Vertex_Triangle_ForVelocity[2] = cuDFNsys::MakeVector2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[2], Coordinate2D_Vec_dev_ptr[EleID - 1].y[2]);
                    /// printf("Error in Particle ID %d, FracID: %d: when the random walker gets into another fracture, the trajectory is not on the same side of the element!\n", i + 1, FracID);
                    /// //printf("Error: the trajectory (particle ID: %d) is acrossing the inlet, but the molecular diffusion is zero\n", i + 1);
                    /// printf("Trajectory of particle ID %d:\n\t[%.40f, %.40f]\n\t[%.40f, %.40f]\n", i + 1,
                    ///        InitPos.x, InitPos.y, TargPos.x, TargPos.y);
                    /// printf("The edge of particle ID %d intersected:\n\t[%.40f, %.40f]\n\t[%.40f, %.40f]\n", i + 1,
                    ///        TheEdge[0].x, TheEdge[0].y, TheEdge[1].x, TheEdge[1].y);
                    /// printf("Intersection of particle ID %d intersected:\n\t[%.40f, %.40f]\n", i + 1,
                    ///        IntersectionOnEdge2D_anotherFrac.x, IntersectionOnEdge2D_anotherFrac.y);
                    /// printf("Element ID %d of particle ID %d:\n\t[%.40f, %.40f]\n\t[%.40f, %.40f]\n\t[%.40f, %.40f]\n", EleID, i + 1,
                    ///        Vertex_Triangle_ForVelocity[0].x, Vertex_Triangle_ForVelocity[0].y, Vertex_Triangle_ForVelocity[1].x, Vertex_Triangle_ForVelocity[1].y,
                    ///        Vertex_Triangle_ForVelocity[2].x, Vertex_Triangle_ForVelocity[2].y);
                    /// cuDFNsys::Vector3<T> InitiPos3D = cuDFNsys::MakeVector3(InitPos.x, InitPos.y, (T)0.0);
                    /// cuDFNsys::Vector3<T> TargPos3D = cuDFNsys::MakeVector3(TargPos.x, TargPos.y, (T)0.0);
                    /// T RK_1[3][3];
                    /// Frac_DEV[FracID].RoationMatrix(RK_1, 23);
                    /// InitiPos3D = cuDFNsys::ProductSquare3Vector3<T>(RK_1, InitiPos3D);
                    /// TargPos3D = cuDFNsys::ProductSquare3Vector3<T>(RK_1, TargPos3D);
                    /// InitiPos3D.x += Frac_DEV[FracID].Center.x;
                    /// InitiPos3D.y += Frac_DEV[FracID].Center.y;
                    /// InitiPos3D.z += Frac_DEV[FracID].Center.z;
                    /// TargPos3D.x += Frac_DEV[FracID].Center.x;
                    /// TargPos3D.y += Frac_DEV[FracID].Center.y;
                    /// TargPos3D.z += Frac_DEV[FracID].Center.z;
                    /// printf("Trajectory3D of ParticleID %d, FracID: %d:\n\t%.40f, %.40f, %.40f\n\t%.40f, %.40f, %.40f\n",
                    ///        i + 1, FracID,
                    ///        InitiPos3D.x,
                    ///        InitiPos3D.y,
                    ///        InitiPos3D.z,
                    ///        TargPos3D.x,
                    ///        TargPos3D.y,
                    ///        TargPos3D.z);
                    /// goto Debug100;
                }
            };

            // P_DEV[i].AccumDisplacement += norm_FinishedTrajec;

            IfUpdateTrajectoryLastStep = false;
            HaveRandomWalkerBeenThroughTrace = true;
            // printf("\nHaveRandomWalkerBeenThroughTrace: true\n");
            normLastTrajectory = pow((TargPos.x - InitPos.x) * (TargPos.x - InitPos.x) + (TargPos.y - InitPos.y) * (TargPos.y - InitPos.y), 0.5);
            continue;
        }
    };
    //-------------------------

    if (i + 1 == 0)
    {
    TimeStepInfo:;
        // Particle_runtime_error_dev_pnt[i] = 1;

        // new point 3D
        cuDFNsys::Vector3<T> TargPos3D;

        if (P_DEV[i].ParticleID >= 0)
            TargPos3D = cuDFNsys::MakeVector3(P_DEV[i].Position2D.x, P_DEV[i].Position2D.y, (T)0.0);
        else
            TargPos3D = cuDFNsys::MakeVector3(TargPos.x, TargPos.y, (T)0.0);

        T RK_1[3][3];
        Frac_DEV[FracID].RoationMatrix(RK_1, 23);
        TargPos3D = cuDFNsys::ProductSquare3Vector3<T>(RK_1, TargPos3D);

        TargPos3D.x += Frac_DEV[FracID].Center.x;
        TargPos3D.y += Frac_DEV[FracID].Center.y;
        TargPos3D.z += Frac_DEV[FracID].Center.z;

        // pre position
        cuDFNsys::Vector3<T> Init_pre_3D = cuDFNsys::MakeVector3(Pre_initPos.x, Pre_initPos.y, (T)0.0);
        Frac_DEV[Pre_fractureIDID].RoationMatrix(RK_1, 23);
        Init_pre_3D = cuDFNsys::ProductSquare3Vector3<T>(RK_1, Init_pre_3D);

        Init_pre_3D.x += Frac_DEV[Pre_fractureIDID].Center.x;
        Init_pre_3D.y += Frac_DEV[Pre_fractureIDID].Center.y;
        Init_pre_3D.z += Frac_DEV[Pre_fractureIDID].Center.z;

        P_DEV[i].AccumDisplacement += sqrt((TargPos3D.x - Init_pre_3D.x) * (TargPos3D.x - Init_pre_3D.x) + (TargPos3D.y - Init_pre_3D.y) * (TargPos3D.y - Init_pre_3D.y) +
                                           (TargPos3D.z - Init_pre_3D.z) * (TargPos3D.z - Init_pre_3D.z));
        //printf("%.30f\n", P_DEV[i].AccumDisplacement);

        return;
    }

    if (i + 1 == 0)
    {
    Debug100:;
        Particle_runtime_error_dev_pnt[i] = 1;
        printf("Particle ID: %d, EleID (from 1): %d, P_DEV[i].Position2D: %.40f, %.40f,\nz1: %.40f\nz2: %.40f\n", i + 1,
               P_DEV[i].ElementID, P_DEV[i].Position2D.x, P_DEV[i].Position2D.y, z1, z2);
        return;
    }
}; // ParticleMovementOneTimeStepGPUKernel
template __global__ void cuDFNsys::ParticleMovementOneTimeStepGPUKernel<double>(unsigned long seed,
                                                                                double delta_T_,
                                                                                double Dispersion_local,
                                                                                cuDFNsys::Particle<double> *P_DEV,
                                                                                cuDFNsys::Fracture<double> *Frac_DEV,
                                                                                cuDFNsys::EdgeToEle *EdgesSharedEle_DEV,
                                                                                cuDFNsys::EleCoor<double> *Coordinate2D_Vec_dev_ptr,
                                                                                //cuDFNsys::NeighborEle *NeighborEleOfOneEle_dev_ptr,
                                                                                uint *EleToFracID_ptr,
                                                                                double *velocity_ptr,
                                                                                uint Dir_flow,
                                                                                double outletcoordinate,
                                                                                int count,
                                                                                int numElements,
                                                                                uint stepNO, uint *Particle_runtime_error_dev_pnt, uint NUMParticlesInTotal,
                                                                                bool If_completeMixing);
template __global__ void cuDFNsys::ParticleMovementOneTimeStepGPUKernel<float>(unsigned long seed,
                                                                               float delta_T_,
                                                                               float Dispersion_local,
                                                                               cuDFNsys::Particle<float> *P_DEV,
                                                                               cuDFNsys::Fracture<float> *Frac_DEV,
                                                                               cuDFNsys::EdgeToEle *EdgesSharedEle_DEV,
                                                                               cuDFNsys::EleCoor<float> *Coordinate2D_Vec_dev_ptr,
                                                                               //cuDFNsys::NeighborEle *NeighborEleOfOneEle_dev_ptr,
                                                                               uint *EleToFracID_ptr,
                                                                               float *velocity_ptr,
                                                                               uint Dir_flow,
                                                                               float outletcoordinate,
                                                                               int count,
                                                                               int numElements,
                                                                               uint stepNO, uint *Particle_runtime_error_dev_pnt, uint NUMParticlesInTotal, bool If_completeMixing);
