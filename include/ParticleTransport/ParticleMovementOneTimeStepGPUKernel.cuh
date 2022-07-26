///////////////////////////////////////////////////////////////////
// NAME:              ParticleMovementOneTimeStepGPUKernel.cuh
//
// PURPOSE:           GPU kernel function to move particles
//                    in one time step
//
// FUNCTIONS/OBJECTS: ParticleMovementOneTimeStepGPUKernel
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../DataTypeSelector/DataTypeSelector.cuh"
#include "../Fractures/Fracture.cuh"
#include "../Geometry/Geometry.cuh"
#include "../GlobalDef/GlobalDef.cuh"
#include "../MHFEM/ReconstructVelocityGrid.cuh"
#include "../MatrixManipulation/MatrixManipulation.cuh"
#include "../RandomFunction/RandomFunction.cuh"
#include "EdgeToEle.cuh"
#include "IdentifyParticleCrossesWhichEdge.cuh"
#include "IfParticleOnBoundOfElement.cuh"
#include "IfParticlePositionInNeighboringElement.cuh"
#include "IfParticlePositionNearOneVertexOfElement.cuh"
#include "Particle.cuh"
#include "ParticleReflection.cuh"
#include "Roate2DPositionTo3D.cuh"
#include "WhichElementToGo.cuh"

namespace cuDFNsys
{
template <typename T>
__global__ void ParticleMovementOneTimeStepGPUKernel(unsigned long seed,
                                                     T delta_T_,
                                                     T Dispersion_local,
                                                     cuDFNsys::Particle<T> *P_DEV,
                                                     cuDFNsys::Fracture<T> *Frac_DEV,
                                                     cuDFNsys::EdgeToEle *EdgesSharedEle_DEV,
                                                     cuDFNsys::EleCoor<T> *Coordinate2D_Vec_dev_ptr,
                                                     cuDFNsys::NeighborEle *NeighborEleOfOneEle_dev_ptr,
                                                     uint *EleToFracID_ptr,
                                                     T *velocity_ptr,
                                                     uint Dir_flow,
                                                     T outletcoordinate,
                                                     int count,
                                                     int numElements,
                                                     uint stepNO)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= count)
        return;

    /// ----------------------- debug -----------------------
    /// ----------------------- debug -----------------------
    /// ----------------------- debug -----------------------
    // P_DEV[i].ElementID = 445;
    // P_DEV[i].Position2D = cuDFNsys::MakeVector2<T>(-0.240469999999999989315213611007, 14.778629999999999711235432187095);
    /// ----------------------- debug -----------------------
    /// ----------------------- debug -----------------------
    /// ----------------------- debug -----------------------

    if (P_DEV[i].IfReachOutletPlane == true)
        return;

    // velocity vector of the grid
    uint EleID = P_DEV[i].ElementID;                                        // from 1
    uint FracID = EleToFracID_ptr[EleID - 1];                               // from 0
    uint3 EdgeNO = make_uint3(EleID * 3 - 3, EleID * 3 - 2, EleID * 3 - 1); // from 0
    cuDFNsys::Vector2<T> InitPos = P_DEV[i].Position2D;
    cuDFNsys::Vector3<T> Veloc_triangle = cuDFNsys::MakeVector3(velocity_ptr[EdgeNO.x],
                                                                velocity_ptr[EdgeNO.y],
                                                                velocity_ptr[EdgeNO.z]);
    cuDFNsys::Vector2<T> Vertex_Triangle_ForVelocity[3];
    Vertex_Triangle_ForVelocity[0] = cuDFNsys::MakeVector2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[0], Coordinate2D_Vec_dev_ptr[EleID - 1].y[0]);
    Vertex_Triangle_ForVelocity[1] = cuDFNsys::MakeVector2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[1], Coordinate2D_Vec_dev_ptr[EleID - 1].y[1]);
    Vertex_Triangle_ForVelocity[2] = cuDFNsys::MakeVector2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[2], Coordinate2D_Vec_dev_ptr[EleID - 1].y[2]);

    cuDFNsys::Vector2<T> Veloc_p = cuDFNsys::ReconstructVelocityGrid<T>(InitPos, Vertex_Triangle_ForVelocity, Veloc_triangle);

    /// ----------------------- debug -----------------------
    /// ----------------------- debug -----------------------
    /// ----------------------- debug -----------------------
    //bool IfTargPosInGrid_yy = cuDFNsys::IfPntInside2DConvexPoly<T>(InitPos, Vertex_Triangle_ForVelocity, 3);
    //if (!IfTargPosInGrid_yy)
    //{
    //    printf("not inside~\n");
    //    return;
    //}
    //
    //printf("velocity 2D: %.30f, %.30f\n", Veloc_p.x, Veloc_p.y);

    /// ----------------------- debug -----------------------
    /// ----------------------- debug -----------------------
    /// ----------------------- debug -----------------------

    // ------------------move the particle
    // ------------------move the particle
    // ------------------move the particle
    curandState state;

    curand_init(seed, i, 0, &state);

    T z1 = curand_normal(&state),
      z2 = curand_normal(&state);

    cuDFNsys::Vector2<T> TargPos;
    TargPos.x = InitPos.x + Veloc_p.x * delta_T_ + (T)1.0f * z1 * sqrt((T)2.0f * Dispersion_local * delta_T_);
    TargPos.y = InitPos.y + Veloc_p.y * delta_T_ + (T)1.0f * z2 * sqrt((T)2.0f * Dispersion_local * delta_T_);

    TargPos.x = round(TargPos.x * 1e5) / 1e5;
    TargPos.y = round(TargPos.y * 1e5) / 1e5;

    /// ----------------------- debug -----------------------
    /// ----------------------- debug -----------------------
    /// ----------------------- debug -----------------------
    // TargPos.x = 0.048109999999999999986677323704;
    // TargPos.y = 15.071410000000000195541360881180;
    /// ----------------------- debug -----------------------
    /// ----------------------- debug -----------------------
    /// ----------------------- debug -----------------------

    //------------record data to debug-------------
    //------------record data to debug-------------
    //------------record data to debug-------------
    cuDFNsys::Vector2<T> designed_particle_trajectory[2] = {InitPos, TargPos};
    uint InitELeID = EleID;
    //---------------------------------------------
    //---------------------------------------------
    //---------------------------------------------

    cuDFNsys::Vector2<T> CrossedGlobalEdge[12][2];
    int CountCrossedGlobalEdge = 0;

    cuDFNsys::Vector2<T> whole_Particle_trajectory[2] = {InitPos, TargPos};

    for (uint Loop_time = 1;; Loop_time++)
    {
        //printf("\n--------------------\nThe looptime: %d\n", Loop_time);

        if (Loop_time >= 10 || CountCrossedGlobalEdge == 11)
        {
            printf("Particle %d, Loop times is too many!\n", i + 1);
            goto Debug100;
        }
        cuDFNsys::Vector2<T> Vertex_Triangle[3];
        Vertex_Triangle[0] = cuDFNsys::MakeVector2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[0], Coordinate2D_Vec_dev_ptr[EleID - 1].y[0]);
        Vertex_Triangle[1] = cuDFNsys::MakeVector2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[1], Coordinate2D_Vec_dev_ptr[EleID - 1].y[1]);
        Vertex_Triangle[2] = cuDFNsys::MakeVector2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[2], Coordinate2D_Vec_dev_ptr[EleID - 1].y[2]);

        // ----------------check if the new position is in the grid----------------
        // ----------------check if the new position is in the grid----------------
        // ----------------check if the new position is in the grid----------------

        bool IfTargPosInGrid = cuDFNsys::IfPntInside2DConvexPoly<T>(TargPos, Vertex_Triangle, 3);
        if (IfTargPosInGrid == true)
        {
            P_DEV[i].ElementID = EleID;
            P_DEV[i].Position2D = TargPos;
            //printf("Particle %d is still in the grid\n", i + 1);
            goto TimeStepInfo;
        }

        // ----------------check if the new position is on one of the bounds of grid?------------------------------
        // ----------------check if the new position is on one of the bounds of grid?------------------------------
        // ----------------check if the new position is on one of the bounds of grid?------------------------------
        int edgelocal = -1;
        bool IfTargPosOnGridBound = cuDFNsys::IfParticleOnBoundOfElement(TargPos,
                                                                         Vertex_Triangle,
                                                                         edgelocal,
                                                                         (T)1e-11);

        if (IfTargPosOnGridBound == true)
        {
            P_DEV[i].ElementID = EleID;
            P_DEV[i].Position2D = TargPos;

            // if this target point reaches the outlet plane?
            cuDFNsys::Vector3<T> Pos3D_ = cuDFNsys::Roate2DPositionTo3D<T>(TargPos,
                                                                           Frac_DEV[EleToFracID_ptr[P_DEV[i].ElementID - 1]]);
            T *tmp_pnt = &(Pos3D_.x);
            if (abs(tmp_pnt[Dir_flow] - outletcoordinate) < 1e-4)
            {
                P_DEV[i].IfReachOutletPlane = true; // reaches outlet plane
                return;
            }
            else
                goto TimeStepInfo;
        };

        ///---------------now, I am sure that the target position is out of the grid very much------------------------
        ///---------------now, I am sure that the target position is out of the grid very much------------------------
        ///---------------now, I am sure that the target position is out of the grid very much------------------------
        ///---------------now, I am sure that the target position is out of the grid very much------------------------
        ///---------------now, I am sure that the target position is out of the grid very much------------------------

        cuDFNsys::Vector2<T> P_trajectory_[2] = {InitPos, TargPos};
        cuDFNsys::Vector3<T> Result_ = cuDFNsys::IdentifyParticleCrossesWhichEdge(P_trajectory_, Vertex_Triangle, (T)1e-7, CrossedGlobalEdge, CountCrossedGlobalEdge, stepNO, i + 1);

        uint GlobalEdgeNO = 0;
        cuDFNsys::Vector2<T> IntersectionOnCrossedEdge;
        if (Result_.z >= 0)
        {
        YK100:;
            IntersectionOnCrossedEdge = cuDFNsys::MakeVector2(Result_.x, Result_.y);
            uint *tmp_k = &(EdgeNO.x);
            GlobalEdgeNO = tmp_k[(uint)Result_.z];
            // printf("\nResult_.z >= 0; Result_.z = %.30f, edgeno: %d, %d, %d\n", Result_.z, EdgeNO.x, EdgeNO.y, EdgeNO.z);

            CrossedGlobalEdge[CountCrossedGlobalEdge][0] = cuDFNsys::MakeVector2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[(uint)Result_.z],
                                                                                 Coordinate2D_Vec_dev_ptr[EleID - 1].y[(uint)Result_.z]);
            CrossedGlobalEdge[CountCrossedGlobalEdge][1] = cuDFNsys::MakeVector2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[((uint)Result_.z + 1) % 3],
                                                                                 Coordinate2D_Vec_dev_ptr[EleID - 1].y[((uint)Result_.z + 1) % 3]);
            CountCrossedGlobalEdge++;
        }
        else
        {
            // cannot find which edge of element does the particle trajectory step over
            // cannot find which edge of element does the particle trajectory step over
            // cannot find which edge of element does the particle trajectory step over

            // this may be due to the precision problem
            // this may be due to the precision problem
            // this may be due to the precision problem

            // let us use the complete particle trajectory!
            // let us use the complete particle trajectory!
            // let us use the complete particle trajectory!

            Result_ = cuDFNsys::IdentifyParticleCrossesWhichEdge(whole_Particle_trajectory,
                                                                 Vertex_Triangle, (T)1e-7,
                                                                 CrossedGlobalEdge, CountCrossedGlobalEdge, stepNO, i + 1);

            if (Result_.z >= 0)
                goto YK100;
            else
            {

                printf("Warning: I cannot find which edge of element %d does the particle trajectory step over! Particle ID: %d\n", EleID, i + 1);
                goto Debug100;
            }
        };

        /// ---------------- if the intersection point is on the inlet plane, and the target point is above the inlet plane? ---------
        /// ---------------- if the intersection point is on the inlet plane, and the target point is above the inlet plane? --------
        /// ---------------- if the intersection point is on the inlet plane, and the target point is above the inlet plane? --------
        /// ---------------- this happens, particularly for the first step when the molecular diffusion is large--------------
        /// ---------------- this happens, particularly for the first step when the molecular diffusion is large--------------
        /// ---------------- this happens, particularly for the first step when the molecular diffusion is large--------------
        cuDFNsys::Vector3<T> IntersectionPnt_3D = cuDFNsys::Roate2DPositionTo3D<T>(IntersectionOnCrossedEdge,
                                                                                   Frac_DEV[EleToFracID_ptr[P_DEV[i].ElementID - 1]]);

        cuDFNsys::Vector3<T> TargetPnt_3D = cuDFNsys::Roate2DPositionTo3D<T>(TargPos,
                                                                             Frac_DEV[EleToFracID_ptr[P_DEV[i].ElementID - 1]]);
        T *tmp4 = &(IntersectionPnt_3D.x);
        if (abs(tmp4[Dir_flow] - ((T)-1.0) * outletcoordinate) < 1e-4)
        {
            tmp4 = &(TargetPnt_3D.x);

            if (tmp4[Dir_flow] < ((T)-1.0) * outletcoordinate && tmp4[Dir_flow] > outletcoordinate)
            {
                // do nothing
            }
            else
            {
                printf("Warning: the particle %d goes above the model or moves along the inlet plane!\n", i + 1);
                // delete this particle in the future
                return;
            }
        }

        /// ----------------------- if the intersection point is on the outlet plane
        /// ----------------------- if the intersection point is on the outlet plane
        /// ----------------------- if the intersection point is on the outlet plane
        tmp4 = &(IntersectionPnt_3D.x);
        if (abs(tmp4[Dir_flow] - outletcoordinate) < 1e-4)
        {
            P_DEV[i].IfReachOutletPlane = true;
            P_DEV[i].ElementID = EleID;
            P_DEV[i].Position2D = IntersectionOnCrossedEdge;
            return;
        }

        uint NumOfElesSharedEdge = EdgesSharedEle_DEV[GlobalEdgeNO].NumSharedEle;
        /// ----------------------- if the intersection point is on a fracture boundary, but nr of sharing edge is > 1
        /// ----------------------- if the intersection point is on a fracture boundary, but nr of sharing edge is > 1
        /// ----------------------- if the intersection point is on a fracture boundary, but nr of sharing edge is > 1
        /// ----------------------- this indicates that the identification of crossed edge is wrong
        /// ----------------------- this indicates that the identification of crossed edge is wrong
        /// ----------------------- this indicates that the identification of crossed edge is wrong
        cuDFNsys::Vector2<T> *Frac2Dvertex = new cuDFNsys::Vector2<T>[Frac_DEV[EleToFracID_ptr[EleID - 1]].NumVertsTruncated];
        Frac_DEV[EleToFracID_ptr[EleID - 1]].Generate2DVerts(Frac2Dvertex, Frac_DEV[EleToFracID_ptr[EleID - 1]].NumVertsTruncated, true);
        int edgeNO = -1;
        bool ifOnFracBound = cuDFNsys::IfPntLiesOnBound2DConvexPolyReturnEdgeNO<T>(IntersectionOnCrossedEdge,
                                                                                   Frac2Dvertex, Frac_DEV[FracID].NumVertsTruncated, (T)1e-4, &edgeNO);
        bool IfTargPosInFrac = cuDFNsys::IfPntInside2DConvexPoly<T>(TargPos, Frac2Dvertex, Frac_DEV[EleToFracID_ptr[EleID - 1]].NumVertsTruncated);

        if (!IfTargPosInFrac && ifOnFracBound && NumOfElesSharedEdge > 1)
        {
            cuDFNsys::Vector2<T> SegEdge[2] = {Frac2Dvertex[edgeNO], Frac2Dvertex[(edgeNO + 1) % Frac_DEV[FracID].NumVertsTruncated]};

            delete[] Frac2Dvertex;
            Frac2Dvertex = NULL;

            InitPos = IntersectionOnCrossedEdge;
            TargPos = cuDFNsys::ParticleReflection<T>(TargPos, SegEdge[0], SegEdge[1]);
            /// because the trajectory is changed, the history crossed edge is clear
            /// because the trajectory is changed, the history crossed edge is clear
            /// because the trajectory is changed, the history crossed edge is clear
            whole_Particle_trajectory[0] = InitPos;
            whole_Particle_trajectory[1] = TargPos;
            CountCrossedGlobalEdge = 0;

            continue;
        }

        delete[] Frac2Dvertex;
        Frac2Dvertex = NULL;

        /// ----------------crossed edge is identified-----------
        /// ----------------crossed edge is identified-----------
        /// ----------------crossed edge is identified-----------
        // printf("\n\nNumOfElesSharedEdge = %d\n\n", NumOfElesSharedEdge);

        if (NumOfElesSharedEdge == 1)
        {
            // we have to check if the edge is the target plane or not
            // we have to check if the edge is the target plane or not
            // we have to check if the edge is the target plane or not
            cuDFNsys::Vector3<T> Test3DIntersectionOnCrossedEdge = cuDFNsys::Roate2DPositionTo3D<T>(IntersectionOnCrossedEdge, Frac_DEV[FracID]);
            T *tmp3 = &(Test3DIntersectionOnCrossedEdge.x);
            if (abs(tmp3[Dir_flow] - outletcoordinate) < 1e-4)
            {
                P_DEV[i].IfReachOutletPlane = true;
                P_DEV[i].ElementID = EleID;
                P_DEV[i].Position2D = IntersectionOnCrossedEdge;
                return;
            };

            // if it is a inlet plane
            // if it is a inlet plane
            // if it is a inlet plane
            if (tmp3[Dir_flow] > ((T)-1.0 * outletcoordinate - 1e-4))
            {
                printf("Warning: the particle %d goes above the model!\n", i + 1);
                // delete this particle in the future
                return;
            }

            // if it is a non-flux plane
            // if it is a non-flux plane
            // if it is a non-flux plane
            if (tmp3[Dir_flow] > outletcoordinate && tmp3[Dir_flow] < ((T)-1.0f * outletcoordinate))
            {
                if (!IfTargPosInFrac)
                {
                    /// the target point is indeed out of fracture
                    /// the target point is indeed out of fracture
                    /// the target point is indeed out of fracture

                    // printf("Particle %d trajectory meets a non flux bound, before reflection, the trajectory is\n\t%.30f, %.30f\n\t%.30f, %.30f\n",
                    //        i + 1, InitPos.x, InitPos.y, TargPos.x, TargPos.y);
                    TargPos = cuDFNsys::ParticleReflection<T>(TargPos, Vertex_Triangle[(uint)Result_.z], Vertex_Triangle[((uint)Result_.z + 1) % 3]);
                    // we do not need to change element ID right now
                    // let us determine element ID in next loop
                    InitPos = IntersectionOnCrossedEdge;
                    // printf("Particle %d trajectory meets a non flux bound, after reflction, the trajectory is\n\t%.30f, %.30f\n\t%.30f, %.30f\n", i + 1, InitPos.x, InitPos.y, TargPos.x, TargPos.y);

                    whole_Particle_trajectory[0] = InitPos;
                    whole_Particle_trajectory[1] = TargPos;
                    // after reflection, the particle may cross again the same edge
                    // after reflection, the particle may cross again the same edge
                    // after reflection, the particle may cross again the same edge
                    // so let's clear the vecter of recording of cross edge
                    // so let's clear the vecter of recording of cross edge
                    // so let's clear the vecter of recording of cross edge
                    //  but keep the non-flux egde!!!
                    //  but keep the non-flux egde!!!
                    //  but keep the non-flux egde!!!
                    CrossedGlobalEdge[0][0] = CrossedGlobalEdge[CountCrossedGlobalEdge - 1][0];
                    CrossedGlobalEdge[0][1] = CrossedGlobalEdge[CountCrossedGlobalEdge - 1][1];

                    CountCrossedGlobalEdge = 1;
                    continue;
                }
                else
                {
                    /// the target point is in the fracture, but the intersection point is just on the fracture bound
                    /// the target point is in the fracture, but the intersection point is just on the fracture bound
                    /// the target point is in the fracture, but the intersection point is just on the fracture bound
                    /// the target point may go out the grid
                    /// the target point may go out the grid
                    /// the target point may go out the grid

                    continue;
                }
            }
        }
        else if (NumOfElesSharedEdge == 2)
        {

            //printf("Particle %d trajectory meets a grid edge, the elementID is %d, the cross global edgeNO is %d, the trajectory is\n\t%.30f, %.30f\n\t%.30f, %.30f\n", i + 1, EleID, EleID * 3 - 3 + (uint)Result_.z, InitPos.x, InitPos.y, TargPos.x, TargPos.y);
            EleID = (EdgesSharedEle_DEV[GlobalEdgeNO].EleID[0] == EleID ? EdgesSharedEle_DEV[GlobalEdgeNO].EleID[1] : EdgesSharedEle_DEV[GlobalEdgeNO].EleID[0]);
            //P_DEV[i].ElementID = EleID;
            InitPos = IntersectionOnCrossedEdge;
            //printf("Particle %d trajectory meets a grid edge, the elementID change to %d, the trajectory is\n\t%.30f, %.30f\n\t%.30f, %.30f\n\tAccumulated crossed edge global NO:\n", i + 1, EleID, InitPos.x, InitPos.y, TargPos.x, TargPos.y);
            FracID = EleToFracID_ptr[EleID - 1];                              // from 0
            EdgeNO = make_uint3(EleID * 3 - 3, EleID * 3 - 2, EleID * 3 - 1); // from 0

            //for (uint y = 0; y < CountCrossedGlobalEdge; ++y)
            //printf("\tParticle %d, crossed global edgeNO:\n\t\t%.30f, %.30f\n\t\t%.30f, %.30f\n", i + 1, CrossedGlobalEdge[y][0].x, CrossedGlobalEdge[y][0].y, CrossedGlobalEdge[y][1].x, CrossedGlobalEdge[y][1].y);
            continue;
        }
        else
        {
            // may go to another fracture
            // may go to another fracture
            // may go to another fracture

            uint PreviousFracID = FracID;

            // determine where to go
            // determine where to go
            // determine where to go
            int IndexLocal = -1; // index of local edge NO the next element;
            // the edge is the one where the intersection lies on

            int newELEID_ = -1;
            int newFracID_ = -1;
            cuDFNsys::WhichElementToGo<T>(EleID,
                                          NumOfElesSharedEdge,
                                          EdgesSharedEle_DEV[GlobalEdgeNO].EleID,
                                          EdgesSharedEle_DEV[GlobalEdgeNO].LocalEdgeNO,
                                          EleToFracID_ptr,
                                          Coordinate2D_Vec_dev_ptr,
                                          velocity_ptr,
                                          curand_uniform(&state),
                                          newELEID_,
                                          newFracID_,
                                          IndexLocal);
            if (newELEID_ == -1 || newFracID_ == -1 || IndexLocal == -1)
            {
                printf("Warning: illeagal indice appear, because the determination of which element to go is failed!\n");
                return;
            }

            FracID = (uint)newFracID_;
            EleID = (uint)newELEID_;
            EdgeNO = make_uint3(EleID * 3 - 3, EleID * 3 - 2, EleID * 3 - 1); // from 0

            if (FracID == PreviousFracID) // did not change fracture ID
            {
                InitPos = IntersectionOnCrossedEdge;
                continue;
            }

            // change initial position
            // change initial position
            // change initial position

            cuDFNsys::Vector3<T> InitiPos3D = cuDFNsys::MakeVector3(IntersectionOnCrossedEdge.x, IntersectionOnCrossedEdge.y, (T)0.0f);

            // printf("before go through intersection, Intersection 2D:\n\t%.30f, %.30f\n", IntersectionOnCrossedEdge.x, IntersectionOnCrossedEdge.y);

            T RK_1[3][3];

            Frac_DEV[PreviousFracID].RoationMatrix(RK_1, 23);
            InitiPos3D = cuDFNsys::ProductSquare3Vector3<T>(RK_1, InitiPos3D);
            InitiPos3D.x += Frac_DEV[PreviousFracID].Center.x;
            InitiPos3D.y += Frac_DEV[PreviousFracID].Center.y;
            InitiPos3D.z += Frac_DEV[PreviousFracID].Center.z;

            // printf("before go through intersection, initial position 3D:\n%.30f, %.30f, %.30f\n\n", InitiPos3D.x, InitiPos3D.y, InitiPos3D.z);

            InitiPos3D.x -= Frac_DEV[FracID].Center.x;
            InitiPos3D.y -= Frac_DEV[FracID].Center.y;
            InitiPos3D.z -= Frac_DEV[FracID].Center.z;
            T RK_2[3][3];
            Frac_DEV[FracID].RoationMatrix(RK_2, 32);
            InitiPos3D = cuDFNsys::ProductSquare3Vector3<T>(RK_2, InitiPos3D);

            InitPos.x = InitiPos3D.x;
            InitPos.y = InitiPos3D.y;

            // change target position
            // change target position
            // change target position
            cuDFNsys::Vector2<T> tmpDis = cuDFNsys::MakeVector2(TargPos.x - IntersectionOnCrossedEdge.x, TargPos.y - IntersectionOnCrossedEdge.y);
            T normDistance = sqrt(tmpDis.x * tmpDis.x + tmpDis.y * tmpDis.y);

            cuDFNsys::Vector3<T> Intersection3D = cuDFNsys::MakeVector3(IntersectionOnCrossedEdge.x, IntersectionOnCrossedEdge.y, (T)0.0f),
                                 Target3D = cuDFNsys::MakeVector3(TargPos.x, TargPos.y, (T)0.0f);

            Intersection3D = cuDFNsys::ProductSquare3Vector3<T>(RK_1, Intersection3D);
            Intersection3D.x += Frac_DEV[PreviousFracID].Center.x;
            Intersection3D.y += Frac_DEV[PreviousFracID].Center.y;
            Intersection3D.z += Frac_DEV[PreviousFracID].Center.z;

            //printf("intersection 3D:\n%.30f, %.30f, %.30f\n\n", Intersection3D.x, Intersection3D.y, Intersection3D.z);

            Target3D = cuDFNsys::ProductSquare3Vector3<T>(RK_1, Target3D);
            Target3D.x += Frac_DEV[PreviousFracID].Center.x;
            Target3D.y += Frac_DEV[PreviousFracID].Center.y;
            Target3D.z += Frac_DEV[PreviousFracID].Center.z;
            // printf("before go through intersection, target position 3D:\n%.30f, %.30f, %.30f\n\n", Target3D.x, Target3D.y, Target3D.z);
            cuDFNsys::Vector3<T> V_ = cuDFNsys::MakeVector3(Target3D.x - Intersection3D.x,
                                                            Target3D.y - Intersection3D.y,
                                                            Target3D.z - Intersection3D.z);
            T norm_V = sqrt(V_.x * V_.x + V_.y * V_.y + V_.z * V_.z);
            V_.x /= norm_V;
            V_.y /= norm_V;
            V_.z /= norm_V;

            ///// (I -nn) * V
            ///// (I -nn) * V
            ///// (I -nn) * V
            cuDFNsys::Vector3<T> V_new = cuDFNsys::ProjectVToPlaneN<T>(V_, Frac_DEV[FracID].NormalVec);

            //// the angle between the V_new and the 3D velocity vector should be < 90?
            //// get the 3D velocity vector first
            //// the angle between the V_new and the 3D velocity vector should be < 90?
            //// get the 3D velocity vector first
            //// the angle between the V_new and the 3D velocity vector should be < 90?
            //// get the 3D velocity vector first

            cuDFNsys::Vector3<T> NewTarPos3D = cuDFNsys::MakeVector3(Intersection3D.x + V_new.x * normDistance - Frac_DEV[FracID].Center.x,
                                                                     Intersection3D.y + V_new.y * normDistance - Frac_DEV[FracID].Center.y,
                                                                     Intersection3D.z + V_new.z * normDistance - Frac_DEV[FracID].Center.z);
            // printf("after go through intersection but no correction, target position 3D:\n%.30f, %.30f, %.30f\n\n", NewTarPos3D.x, NewTarPos3D.y, NewTarPos3D.z);
            NewTarPos3D = cuDFNsys::ProductSquare3Vector3<T>(RK_2, NewTarPos3D);
            cuDFNsys::Vector2<T> newTagPos2D = cuDFNsys::MakeVector2(NewTarPos3D.x, NewTarPos3D.y);

            // test if the particle goes in an opposite way?
            // test if the particle goes in an opposite way?
            // test if the particle goes in an opposite way?
            cuDFNsys::Vector2<T> EdgeSeg[2];
            EdgeSeg[0].x = Coordinate2D_Vec_dev_ptr[EleID - 1].x[IndexLocal], EdgeSeg[0].y = Coordinate2D_Vec_dev_ptr[EleID - 1].y[IndexLocal];
            EdgeSeg[1].x = Coordinate2D_Vec_dev_ptr[EleID - 1].x[(IndexLocal + 1) % 3], EdgeSeg[1].y = Coordinate2D_Vec_dev_ptr[EleID - 1].y[(IndexLocal + 1) % 3];

            // because FracID changed, the vector of CrossedGlobalEdge should be cleared
            // because FracID changed, the vector of CrossedGlobalEdge should be cleared
            // because FracID changed, the vector of CrossedGlobalEdge should be cleared
            CrossedGlobalEdge[0][0] = EdgeSeg[0];
            CrossedGlobalEdge[0][1] = EdgeSeg[1];

            CountCrossedGlobalEdge = 1;

            cuDFNsys::Vector2<T> CenterGrid;
            CenterGrid.x = (Coordinate2D_Vec_dev_ptr[EleID - 1].x[0] + Coordinate2D_Vec_dev_ptr[EleID - 1].x[1] + Coordinate2D_Vec_dev_ptr[EleID - 1].x[2]) / 3.0f;
            CenterGrid.y = (Coordinate2D_Vec_dev_ptr[EleID - 1].y[0] + Coordinate2D_Vec_dev_ptr[EleID - 1].y[1] + Coordinate2D_Vec_dev_ptr[EleID - 1].y[2]) / 3.0f;

            // printf("no correction, ele to %d, 2D trajectory:\n\t%.30f, %.30f\n\t%.30f, %.30f\n", EleID, InitPos.x, InitPos.y, newTagPos2D.x, newTagPos2D.y);
            uint O1 = cuDFNsys::OrientationThree2DPnts<T>(EdgeSeg[0], EdgeSeg[1], newTagPos2D, (T)1e-7);
            uint O2 = cuDFNsys::OrientationThree2DPnts<T>(EdgeSeg[0], EdgeSeg[1], CenterGrid, (T)1e-7);

            // {
            //     printf("Particle %d goes to another fracture, from fractureID %d to %d, length of remainning trajectory = %.30f\n", i + 1, PreviousFracID, newFracID_, normDistance);
            //     printf("Particle %d, intersection on the crossed edge in 3D: %.30f, %.30f, %.30f\n", i + 1, Intersection3D.x, Intersection3D.y, Intersection3D.z);
            //     printf("Particle %d, the previous velocity direction vector: %.30f, %.30f, %.30f in frac %d\n", i + 1, V_.x, V_.y, V_.z, PreviousFracID);
            //     printf("Particle %d, the new velocity direction vector: %.30f, %.30f, %.30f in frac %d\n", i + 1, V_new.x, V_new.y, V_new.z, newFracID_);
            //     printf("Particle %d, the length of remainning trajectory: %.30f\n", i + 1, normDistance);
            // }

            if (O1 == 0)
            {
                printf("Warning: particle trajectory is along an intersection!\n");
                return;
            }

            if (O1 != 0 && O1 == O2)
            {
                TargPos = newTagPos2D;
                whole_Particle_trajectory[0] = InitPos, whole_Particle_trajectory[1] = TargPos;
                //return;

                // cuDFNsys::Vector2<T> Vertex_Triangle_PPP[3];
                // Vertex_Triangle_PPP[0] = cuDFNsys::MakeVector2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[0], Coordinate2D_Vec_dev_ptr[EleID - 1].y[0]);
                // Vertex_Triangle_PPP[1] = cuDFNsys::MakeVector2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[1], Coordinate2D_Vec_dev_ptr[EleID - 1].y[1]);
                // Vertex_Triangle_PPP[2] = cuDFNsys::MakeVector2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[2], Coordinate2D_Vec_dev_ptr[EleID - 1].y[2]);
                // printf("O1 == O2, ParticleID %d, eleID: %d;\nPntTrajectory:\n%f, %f\n%f, %f\nTri:\n%f, %f\n%f, %f\n%f, %f\nIntersection3D:\n%f, %f, %f\nV_new:\n%f, %f, %f\nNewTarPos3D:\n%f, %f, %f\n",
                //        i + 1, EleID,
                //        InitPos.x, InitPos.y, TargPos.x, TargPos.y,
                //        Vertex_Triangle_PPP[0].x, Vertex_Triangle_PPP[0].y,
                //        Vertex_Triangle_PPP[1].x, Vertex_Triangle_PPP[1].y,
                //        Vertex_Triangle_PPP[2].x, Vertex_Triangle_PPP[2].y,
                //        Intersection3D.x, Intersection3D.y, Intersection3D.z,
                //        V_new.x, V_new.y, V_new.z,
                //        NewTarPos3D.x, NewTarPos3D.y, NewTarPos3D.z);

                continue;
            }
            else
            {
                NewTarPos3D = cuDFNsys::MakeVector3(Intersection3D.x - V_new.x * normDistance - Frac_DEV[FracID].Center.x,
                                                    Intersection3D.y - V_new.y * normDistance - Frac_DEV[FracID].Center.y,
                                                    Intersection3D.z - V_new.z * normDistance - Frac_DEV[FracID].Center.z);
                //printf("after correction, target position 3D:\n%.30f, %.30f, %.30f\n\n", NewTarPos3D.x, NewTarPos3D.y, NewTarPos3D.z);
                NewTarPos3D = cuDFNsys::ProductSquare3Vector3<T>(RK_2, NewTarPos3D);
                newTagPos2D = cuDFNsys::MakeVector2(NewTarPos3D.x, NewTarPos3D.y);
                TargPos = newTagPos2D;
                whole_Particle_trajectory[0] = InitPos, whole_Particle_trajectory[1] = TargPos;
                //return;

                // cuDFNsys::Vector2<T> Vertex_Triangle_PPP[3];
                // Vertex_Triangle_PPP[0] = cuDFNsys::MakeVector2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[0], Coordinate2D_Vec_dev_ptr[EleID - 1].y[0]);
                // Vertex_Triangle_PPP[1] = cuDFNsys::MakeVector2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[1], Coordinate2D_Vec_dev_ptr[EleID - 1].y[1]);
                // Vertex_Triangle_PPP[2] = cuDFNsys::MakeVector2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[2], Coordinate2D_Vec_dev_ptr[EleID - 1].y[2]);
                //printf("O1 != O2, ParticleID %d, eleID: %d;\nPntTrajectory:\n%f, %f\n%f, %f\nTri:\n%f, %f\n%f, %f\n%f, %f\nIntersection3D:\n%f, %f, %f\nV_new:\n%f, %f, %f\nNewTarPos3D:\n%f, %f, %f\n",
                //       i + 1, EleID,
                //       InitPos.x, InitPos.y, TargPos.x, TargPos.y,
                //       Vertex_Triangle_PPP[0].x, Vertex_Triangle_PPP[0].y,
                //       Vertex_Triangle_PPP[1].x, Vertex_Triangle_PPP[1].y,
                //       Vertex_Triangle_PPP[2].x, Vertex_Triangle_PPP[2].y,
                //       Intersection3D.x, Intersection3D.y, Intersection3D.z,
                //       V_new.x, V_new.y, V_new.z,
                //       NewTarPos3D.x, NewTarPos3D.y, NewTarPos3D.z);
                continue;
            }
        };
    };

    if (i + 1 == -1)
    {
    Debug100:;
        {
            printf("Warning: dropped time step!!! ParticleID: %d\n%%Trajectory:\n\t%.30f, %.30f\n\t%.30f, %.30f\nPre-element ID %d, fracID: %d, coordinates:\n\t%.30f, %.30f\n\t%.30f, %.30f\n\t%.30f, %.30f\n", i + 1,
                   designed_particle_trajectory[0].x, designed_particle_trajectory[0].y,
                   designed_particle_trajectory[1].x, designed_particle_trajectory[1].y,
                   InitELeID, EleToFracID_ptr[InitELeID - 1],
                   Vertex_Triangle_ForVelocity[0].x, Vertex_Triangle_ForVelocity[0].y,
                   Vertex_Triangle_ForVelocity[1].x, Vertex_Triangle_ForVelocity[1].y,
                   Vertex_Triangle_ForVelocity[2].x, Vertex_Triangle_ForVelocity[2].y);

            for (uint k = 0; k < NeighborEleOfOneEle_dev_ptr[InitELeID - 1].NumNeighborEle; ++k)
            {
                uint ele2 = NeighborEleOfOneEle_dev_ptr[InitELeID - 1].EleID[k];

                cuDFNsys::Vector2<T> Vertex_Triangle_PPP[3];
                Vertex_Triangle_PPP[0] = cuDFNsys::MakeVector2(Coordinate2D_Vec_dev_ptr[ele2 - 1].x[0], Coordinate2D_Vec_dev_ptr[ele2 - 1].y[0]);
                Vertex_Triangle_PPP[1] = cuDFNsys::MakeVector2(Coordinate2D_Vec_dev_ptr[ele2 - 1].x[1], Coordinate2D_Vec_dev_ptr[ele2 - 1].y[1]);
                Vertex_Triangle_PPP[2] = cuDFNsys::MakeVector2(Coordinate2D_Vec_dev_ptr[ele2 - 1].x[2], Coordinate2D_Vec_dev_ptr[ele2 - 1].y[2]);
                printf("%% ParticleID: %d, neighboring elementID: %d, fracID: %d:\n\t%.30f, %.30f\n\t%.30f, %.30f\n\t%.30f, %.30f\n",
                       i + 1, ele2, EleToFracID_ptr[ele2 - 1],
                       Vertex_Triangle_PPP[0].x, Vertex_Triangle_PPP[0].y,
                       Vertex_Triangle_PPP[1].x, Vertex_Triangle_PPP[1].y,
                       Vertex_Triangle_PPP[2].x, Vertex_Triangle_PPP[2].y);
            }
        };
        return;
    }

    if (i + 1 == -1)
    {
    TimeStepInfo:;
        return;
        // {
        //     cuDFNsys::Vector2<T> Vertex_Triangle_PPP[3];
        //     Vertex_Triangle_PPP[0] = cuDFNsys::MakeVector2(Coordinate2D_Vec_dev_ptr[P_DEV[i].ElementID - 1].x[0], Coordinate2D_Vec_dev_ptr[P_DEV[i].ElementID - 1].y[0]);
        //     Vertex_Triangle_PPP[1] = cuDFNsys::MakeVector2(Coordinate2D_Vec_dev_ptr[P_DEV[i].ElementID - 1].x[1], Coordinate2D_Vec_dev_ptr[P_DEV[i].ElementID - 1].y[1]);
        //     Vertex_Triangle_PPP[2] = cuDFNsys::MakeVector2(Coordinate2D_Vec_dev_ptr[P_DEV[i].ElementID - 1].x[2], Coordinate2D_Vec_dev_ptr[P_DEV[i].ElementID - 1].y[2]);
        //     printf("Time step infomation. ParticleID: %d\nTrajectory:\n\t%.30f, %.30f\n\t%.30f, %.30f\nPre-element ID %d, coordinates:\n\t%.30f, %.30f\n\t%.30f, %.30f\n\t%.30f, %.30f\nParticleID: %d, Now, element ID: %d, coordinates:\n\t%.30f, %.30f\n\t%.30f, %.30f\n\t%.30f, %.30f\n",
        //            i + 1,
        //            designed_particle_trajectory[0].x, designed_particle_trajectory[0].y,
        //            designed_particle_trajectory[1].x, designed_particle_trajectory[1].y,
        //            InitELeID,
        //            Vertex_Triangle_ForVelocity[0].x, Vertex_Triangle_ForVelocity[0].y,
        //            Vertex_Triangle_ForVelocity[1].x, Vertex_Triangle_ForVelocity[1].y,
        //            Vertex_Triangle_ForVelocity[2].x, Vertex_Triangle_ForVelocity[2].y,
        //            i + 1, P_DEV[i].ElementID,
        //            Vertex_Triangle_PPP[0].x, Vertex_Triangle_PPP[0].y,
        //            Vertex_Triangle_PPP[1].x, Vertex_Triangle_PPP[1].y,
        //            Vertex_Triangle_PPP[2].x, Vertex_Triangle_PPP[2].y);
        //     for (uint k = 0; k < NeighborEleOfOneEle_dev_ptr[InitELeID - 1].NumNeighborEle; ++k)
        //     {
        //         uint ele2 = NeighborEleOfOneEle_dev_ptr[InitELeID - 1].EleID[k];
        //         cuDFNsys::Vector2<T> Vertex_Triangle_PPP[3];
        //         Vertex_Triangle_PPP[0] = cuDFNsys::MakeVector2(Coordinate2D_Vec_dev_ptr[ele2 - 1].x[0], Coordinate2D_Vec_dev_ptr[ele2 - 1].y[0]);
        //         Vertex_Triangle_PPP[1] = cuDFNsys::MakeVector2(Coordinate2D_Vec_dev_ptr[ele2 - 1].x[1], Coordinate2D_Vec_dev_ptr[ele2 - 1].y[1]);
        //         Vertex_Triangle_PPP[2] = cuDFNsys::MakeVector2(Coordinate2D_Vec_dev_ptr[ele2 - 1].x[2], Coordinate2D_Vec_dev_ptr[ele2 - 1].y[2]);
        //         printf("ParticleID: %d, neighboring elementID: %d:\n\t%.30f, %.30f\n\t%.30f, %.30f\n\t%.30f, %.30f\n",
        //                i + 1, ele2,
        //                Vertex_Triangle_PPP[0].x, Vertex_Triangle_PPP[0].y,
        //                Vertex_Triangle_PPP[1].x, Vertex_Triangle_PPP[1].y,
        //                Vertex_Triangle_PPP[2].x, Vertex_Triangle_PPP[2].y);
        //     }
        // }
        // return;
    }
};
template __global__ void ParticleMovementOneTimeStepGPUKernel<double>(unsigned long seed,
                                                                      double delta_T_,
                                                                      double Dispersion_local,
                                                                      cuDFNsys::Particle<double> *P_DEV,
                                                                      cuDFNsys::Fracture<double> *Frac_DEV,
                                                                      cuDFNsys::EdgeToEle *EdgesSharedEle_DEV,
                                                                      cuDFNsys::EleCoor<double> *Coordinate2D_Vec_dev_ptr,
                                                                      cuDFNsys::NeighborEle *NeighborEleOfOneEle_dev_ptr,
                                                                      uint *EleToFracID_ptr,
                                                                      double *velocity_ptr,
                                                                      uint Dir_flow,
                                                                      double outletcoordinate,
                                                                      int count,
                                                                      int numElements,
                                                                      uint stepNO);
template __global__ void ParticleMovementOneTimeStepGPUKernel<float>(unsigned long seed,
                                                                     float delta_T_,
                                                                     float Dispersion_local,
                                                                     cuDFNsys::Particle<float> *P_DEV,
                                                                     cuDFNsys::Fracture<float> *Frac_DEV,
                                                                     cuDFNsys::EdgeToEle *EdgesSharedEle_DEV,
                                                                     cuDFNsys::EleCoor<float> *Coordinate2D_Vec_dev_ptr,
                                                                     cuDFNsys::NeighborEle *NeighborEleOfOneEle_dev_ptr,
                                                                     uint *EleToFracID_ptr,
                                                                     float *velocity_ptr,
                                                                     uint Dir_flow,
                                                                     float outletcoordinate,
                                                                     int count,
                                                                     int numElements,
                                                                     uint stepNO);
}; // namespace cuDFNsys
