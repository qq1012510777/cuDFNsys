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
#include "../Geometry/2D/DistancePnt2DSeg.cuh"
#include "../Geometry/2D/IfPntInside2DConvexPoly.cuh"
#include "../Geometry/2D/IfPntLiesOnBound2DConvexPoly.cuh"
#include "../Geometry/2D/IntersectionTwo2DSegs.cuh"
#include "../MHFEM/MHFEM.cuh"
#include "../MHFEM/ReconstructVelocityGrid.cuh"
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
__global__ void ParticleMovementOneTimeStepGPUKernel(unsigned long seed,
                                                     float delta_T_,
                                                     float Dispersion_local,
                                                     cuDFNsys::Particle *P_DEV,
                                                     cuDFNsys::Fracture *Frac_DEV,
                                                     cuDFNsys::EdgeToEle *EdgesSharedEle_DEV,
                                                     cuDFNsys::EleCoor *Coordinate2D_Vec_dev_ptr,
                                                     cuDFNsys::NeighborEle *NeighborEleOfOneEle_dev_ptr,
                                                     uint *EleToFracID_ptr,
                                                     float *velocity_ptr,
                                                     uint Dir_flow,
                                                     float outletcoordinate,
                                                     int count,
                                                     int numElements,
                                                     uint stepNO)
{
    //Restart_RW:;
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= count)
        return;

    if (P_DEV[i].IfReachOutletPlane == true)
        return;

    // velocity vector of the grid
    uint EleID = P_DEV[i].ElementID;                                        // from 1
    uint FracID = EleToFracID_ptr[EleID - 1];                               // from 0
    uint3 EdgeNO = make_uint3(EleID * 3 - 3, EleID * 3 - 2, EleID * 3 - 1); // from 0
    float2 InitPos = P_DEV[i].Position2D;
    float3 Veloc_triangle = make_float3(velocity_ptr[EdgeNO.x],
                                        velocity_ptr[EdgeNO.y],
                                        velocity_ptr[EdgeNO.z]);
    float2 Vertex_Triangle_ForVelocity[3];
    Vertex_Triangle_ForVelocity[0] = make_float2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[0], Coordinate2D_Vec_dev_ptr[EleID - 1].y[0]);
    Vertex_Triangle_ForVelocity[1] = make_float2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[1], Coordinate2D_Vec_dev_ptr[EleID - 1].y[1]);
    Vertex_Triangle_ForVelocity[2] = make_float2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[2], Coordinate2D_Vec_dev_ptr[EleID - 1].y[2]);

    float2 Veloc_p = cuDFNsys::ReconstructVelocityGrid(InitPos, Vertex_Triangle_ForVelocity, Veloc_triangle);

    // ------------------move the particle
    // ------------------move the particle
    // ------------------move the particle
    curandState state;

    curand_init(seed, i, 0, &state);

    float z1 = curand_normal(&state),
          z2 = curand_normal(&state);

    float2 TargPos;
    TargPos.x = InitPos.x + Veloc_p.x * delta_T_ + 1.0f * z1 * sqrt(2.0f * Dispersion_local * delta_T_);
    TargPos.y = InitPos.y + Veloc_p.y * delta_T_ + 1.0f * z2 * sqrt(2.0f * Dispersion_local * delta_T_);

    //------------record data to debug-------------
    //------------record data to debug-------------
    //------------record data to debug-------------
    float2 designed_particle_trajectory[2] = {InitPos, TargPos};
    uint InitELeID = EleID;

    float2 CrossedGlobalEdge[10][2];
    int CountCrossedGlobalEdge = 0;

    float2 whole_Particle_trajectory[2] = {InitPos, TargPos};

    for (uint Loop_time = 1;; Loop_time++)
    {
        if (Loop_time >= 8 || CountCrossedGlobalEdge == 9)
        {
            //printf("Particle %d, Loop times is too many!\n", i + 1);
            goto Debug100;
        }
        float2 Vertex_Triangle[3];
        Vertex_Triangle[0] = make_float2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[0], Coordinate2D_Vec_dev_ptr[EleID - 1].y[0]);
        Vertex_Triangle[1] = make_float2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[1], Coordinate2D_Vec_dev_ptr[EleID - 1].y[1]);
        Vertex_Triangle[2] = make_float2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[2], Coordinate2D_Vec_dev_ptr[EleID - 1].y[2]);

        // ----------------check if the new position is in the grid----------------
        // ----------------check if the new position is in the grid----------------
        // ----------------check if the new position is in the grid----------------

        bool IfTargPosInGrid = cuDFNsys::IfPntInside2DConvexPoly(TargPos, Vertex_Triangle, 3);
        if (IfTargPosInGrid == true)
        {
            P_DEV[i].ElementID = EleID;
            P_DEV[i].Position2D = TargPos;
            //printf("Particle %d is still in the grid\n", i + 1);
            goto TimeStepInfo;
        }

        // ----------------check if the new position is very very close to one vertex of the grid but out of grid?----------------
        // ----------------check if the new position is very very close to one vertex of the grid but out of grid?----------------
        // ----------------check if the new position is very very close to one vertex of the grid but out of grid?----------------
        bool IfTargPosNearVertex = cuDFNsys::IfParticlePositionNearOneVertexOfElement(TargPos, Vertex_Triangle, 1e-7);

        if (IfTargPosNearVertex == true)
        {
            // the position is out of the grid but very very near one vertex of the grid
            // the position is out of the grid but very very near one vertex of the grid
            // the position is out of the grid but very very near one vertex of the grid

            bool GRWQ = cuDFNsys::IfParticlePositionInNeighboringElement(TargPos, EleID, EleToFracID_ptr, FracID,
                                                                         NeighborEleOfOneEle_dev_ptr[EleID - 1], Coordinate2D_Vec_dev_ptr, 1e-11);

            if (GRWQ)
            {
                P_DEV[i].ElementID = EleID;
                P_DEV[i].Position2D = TargPos;
                //printf("Particle %d new position is near a vertex of and out of the grid, so I find it is in the neighboring element is %d. ", i + 1, EleID);
                goto TimeStepInfo;
            }
            else
            {
                printf("Warning: I cannot find which neighboring element does the target position fall into! The target position is out of the grid but very near one of the vertexes of the grid.\n");
                goto Debug100;
            }
        };

        // ----------------check if the new position is on one of the bounds of grid?------------------------------
        // ----------------check if the new position is on one of the bounds of grid?------------------------------
        // ----------------check if the new position is on one of the bounds of grid?------------------------------
        int edgelocal = -1;
        bool IfTargPosOnGridBound = cuDFNsys::IfParticleOnBoundOfElement(TargPos,
                                                                         Vertex_Triangle,
                                                                         edgelocal,
                                                                         1e-7);

        if (IfTargPosOnGridBound == true)
        {
            // since the target position is not inside the element/grid, the target position must outside the element or lie on the bound, which is just very very very close to bound
            // since the target position is not inside the element/grid, the target position must outside the element or lie on the bound, which is just very very very close to bound
            // since the target position is not inside the element/grid, the target position must outside the element or lie on the bound, which is just very very very close to bound

            // determine which element to go now!!!
            // determine which element to go now!!!
            // determine which element to go now!!!
            //printf("Target position is extramely close to one bound of grid\n");

            // but how about if the target position is really on one bound of grid???
            // but how about if the target position is really on one bound of grid???
            // but how about if the target position is really on one bound of grid???

            float2 EdgeH[2];
            EdgeH[0] = Vertex_Triangle[edgelocal];
            EdgeH[1] = Vertex_Triangle[(edgelocal + 1) % 3];

            float Distance_Seg_Tar = cuDFNsys::DistancePnt2DSeg(TargPos, EdgeH);

            if (abs(Distance_Seg_Tar) < 1e-13)
            {
                // this error is very very small!!!
                P_DEV[i].ElementID = EleID;
                P_DEV[i].Position2D = TargPos;
                goto TimeStepInfo;
            }
            // else go to next piece of code
            // else go to next piece of code
            // else go to next piece of code
        };

        ///---------------now, I am sure that the target position is out of the grid very much------------------------
        ///---------------now, I am sure that the target position is out of the grid very much------------------------
        ///---------------now, I am sure that the target position is out of the grid very much------------------------
        ///---------------now, I am sure that the target position is out of the grid very much------------------------
        ///---------------now, I am sure that the target position is out of the grid very much------------------------

        if (!IfTargPosOnGridBound && !IfTargPosNearVertex)
        {
            ///-------------if the initial position is just near one of the vertexes of the grid while the target position is out of the grid--------------------
            ///-------------if the initial position is just near one of the vertexes of the grid while the target position is out of the grid--------------------
            ///-------------if the initial position is just near one of the vertexes of the grid while the target position is out of the grid--------------------

            bool IfInitPosNearVertex = cuDFNsys::IfParticlePositionNearOneVertexOfElement(InitPos, Vertex_Triangle, 1e-4);

            if (IfInitPosNearVertex)
            {
                //------------test if the initial position is on target plane or not---------------
                //------------test if the initial position is on target plane or not---------------
                //------------test if the initial position is on target plane or not---------------
                float3 Test3DInitPos = cuDFNsys::Roate2DPositionTo3D(InitPos, Frac_DEV[FracID]);
                float *Tmp = &(Test3DInitPos.x);
                if (abs(Tmp[Dir_flow] - outletcoordinate) < 1e-3)
                {
                    // actually the initial position is already on the target plane, we can discontinue the simulation
                    P_DEV[i].IfReachOutletPlane = true;
                    return;
                }
                /// now, the initial position could be on non flux boundary or inside the frac...
                /// now, the initial position could be on non flux boundary or inside the frac...
                /// now, the initial position could be on non flux boundary or inside the frac...

                // -------------test a special case: if both the vertex and the initial position are on the non-flux edge!!!-----------
                // -------------test a special case: if both the vertex and the initial position are on the non-flux edge!!!-----------
                // -------------test a special case: if both the vertex and the initial position are on the non-flux edge!!!-----------

                // namely, on one of bounds of the fracture
                int edgeNO = -1;
                float2 *Frac2Dvertex = new float2[Frac_DEV[FracID].NumVertsTruncated];
                Frac_DEV[FracID].Generate2DVerts(Frac2Dvertex, Frac_DEV[FracID].NumVertsTruncated, true);
                bool IfOnFracBound = cuDFNsys::IfPntLiesOnBound2DConvexPolyReturnEdgeNO(InitPos, Frac2Dvertex, Frac_DEV[FracID].NumVertsTruncated, 1e-4, &edgeNO);

                if (IfOnFracBound)
                {
                    //---------------------------if the intial position is on bound and the target position goes above the model------------------------
                    //---------------------------if the intial position is on bound and the target position goes above the model------------------------
                    //---------------------------if the intial position is on bound and the target position goes above the model------------------------
                    float3 Test3DTargPos_1 = cuDFNsys::Roate2DPositionTo3D(TargPos, Frac_DEV[FracID]);
                    float *tmp2 = &(Test3DTargPos_1.x);
                    if (tmp2[Dir_flow] > -outletcoordinate)
                    {
                        //
                        //printf("Warning: the particle %d goes above the model! Initial position is near a vertex\n", i + 1);
                        delete[] Frac2Dvertex;
                        Frac2Dvertex = NULL;
                        return;
                    };

                    //---------------------------if the intial position is on bound and the target position is out of the frac------------------------
                    //---------------------------if the intial position is on bound and the target position is out of the frac------------------------
                    //---------------------------if the intial position is on bound and the target position is out of the frac------------------------
                    bool IfTargPosOnFracbound = cuDFNsys::IfPntLiesOnBound2DConvexPoly(TargPos, Frac2Dvertex, Frac_DEV[FracID].NumVertsTruncated, 1e-4);
                    bool IfTargPosInFrac = cuDFNsys::IfPntInside2DConvexPoly(TargPos, Frac2Dvertex, Frac_DEV[FracID].NumVertsTruncated);

                    float2 SegEdge[2] = {Frac2Dvertex[edgeNO], Frac2Dvertex[(edgeNO + 1) % Frac_DEV[FracID].NumVertsTruncated]};
                    delete[] Frac2Dvertex;
                    Frac2Dvertex = NULL;

                    if (IfTargPosOnFracbound == false && IfTargPosInFrac == false)
                    {
                        // do reflection about edge NO (edgeNO)
                        // do reflection about edge NO (edgeNO)
                        // do reflection about edge NO (edgeNO)
                        TargPos = cuDFNsys::ParticleReflection(TargPos, SegEdge[0], SegEdge[1]);

                        whole_Particle_trajectory[0] = InitPos;
                        whole_Particle_trajectory[1] = TargPos;

                        // I assume that the particle trajectory is small enough, so that the target position must fall into one neighboring element
                        // I assume that the particle trajectory is small enough, so that the target position must fall into one neighboring element
                        // I assume that the particle trajectory is small enough, so that the target position must fall into one neighboring element

                        // may go back to the same grid!
                        // may go back to the same grid!
                        // may go back to the same grid!
                        bool IfInGrid = cuDFNsys::IfPntInside2DConvexPoly(TargPos, Vertex_Triangle, 3);

                        if (IfInGrid == true)
                        {
                            P_DEV[i].ElementID = EleID;
                            P_DEV[i].Position2D = TargPos;
                            //printf("Particle %d. The initial position is near a vertex and also on the non-flux bound; after reflection, the new position is still in element %d\n", i + 1, EleID);
                            goto TimeStepInfo;
                        }

                        //neighboring element
                        //neighboring element
                        //neighboring element
                        bool IfTargPosInNeighboringElement = cuDFNsys::IfParticlePositionInNeighboringElement(TargPos, EleID, EleToFracID_ptr,
                                                                                                              FracID, NeighborEleOfOneEle_dev_ptr[EleID - 1],
                                                                                                              Coordinate2D_Vec_dev_ptr, 1e-11);

                        if (IfTargPosInNeighboringElement)
                        {
                            P_DEV[i].ElementID = EleID;
                            P_DEV[i].Position2D = TargPos;
                            //printf("Particle %d. The initial position is near a vertex and also on the non-flux bound; after reflection, the new position is ina neighboring element %d\n", i + 1, EleID);
                            goto TimeStepInfo;
                        }
                        else
                        {
                            /// may be the target position reaches target plane
                            /// may be the target position reaches target plane
                            /// may be the target position reaches target plane
                            float3 Test3DTargPos_2 = cuDFNsys::Roate2DPositionTo3D(TargPos, Frac_DEV[FracID]);
                            float *tmp3 = &(Test3DTargPos_2.x);
                            if (tmp3[Dir_flow] <= outletcoordinate)
                            {
                                P_DEV[i].IfReachOutletPlane = true;
                                P_DEV[i].ElementID = EleID;
                                P_DEV[i].Position2D = TargPos;
                                return;
                            };

                            printf("Warning: after reflection, I cannot find which neighboring element does the target positon fall into; particle ID is %d. Initial position is near a vertex\n", i + 1);
                            goto Debug100;
                        }
                    };
                }
                //---------------- if the initial position is near a vertex but not near the non-flux bound-------------
                //---------------- if the initial position is near a vertex but not near the non-flux bound-------------
                //---------------- if the initial position is near a vertex but not near the non-flux bound-------------
                // --------------- general case
                // --------------- general case
                // --------------- general case
                bool IfTargPosInNeighboringElement = cuDFNsys::IfParticlePositionInNeighboringElement(TargPos, EleID,
                                                                                                      EleToFracID_ptr,
                                                                                                      FracID, NeighborEleOfOneEle_dev_ptr[EleID - 1],
                                                                                                      Coordinate2D_Vec_dev_ptr, 1e-4);
                if (IfTargPosInNeighboringElement)
                {
                    P_DEV[i].ElementID = EleID;
                    P_DEV[i].Position2D = TargPos;
                    //printf("Particle %d, the initial position is near a vertex, and in general case (not near a bound), I found the target position is in neighboring element %d\n", i + 1, EleID);
                    goto TimeStepInfo;
                }
                else
                {
                    /// may be the target position reaches target plane
                    /// may be the target position reaches target plane
                    /// may be the target position reaches target plane
                    float3 Test3DTargPos_3 = cuDFNsys::Roate2DPositionTo3D(TargPos, Frac_DEV[FracID]);
                    float *tmp3 = &(Test3DTargPos_3.x);
                    if (tmp3[Dir_flow] < (outletcoordinate + 1e-4))
                    {
                        P_DEV[i].IfReachOutletPlane = true;
                        P_DEV[i].ElementID = EleID;
                        P_DEV[i].Position2D = TargPos;
                        return;
                    };

                    printf("Warning: I cannot find which neighboring element does the target positon fall into; particle ID is %d. Initial position is near a vertex\n", i + 1);
                    goto Debug100;
                }
            }
        }

        ///------------- initial position is not near a vertex------------
        ///------------- initial position is not near a vertex------------
        ///------------- initial position is not near a vertex------------

        float2 P_trajectory_[2] = {InitPos, TargPos};
        float3 Result_ = cuDFNsys::IdentifyParticleCrossesWhichEdge(P_trajectory_, Vertex_Triangle, 1e-7, CrossedGlobalEdge, CountCrossedGlobalEdge, stepNO, i + 1);

        uint GlobalEdgeNO = 0;
        float2 IntersectionOnCrossedEdge;
        if (Result_.z >= 0)
        {
        YK100:;
            IntersectionOnCrossedEdge = make_float2(Result_.x, Result_.y);
            uint *tmp_k = &(EdgeNO.x);
            GlobalEdgeNO = tmp_k[(uint)Result_.z];
            //printf("Result_.z >= 0; Result_.z = %.30f, edgeno: %d, %d, %d\n", Result_.z, EdgeNO.x, EdgeNO.y, EdgeNO.z);

            CrossedGlobalEdge[CountCrossedGlobalEdge][0] = make_float2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[(uint)Result_.z],
                                                                       Coordinate2D_Vec_dev_ptr[EleID - 1].y[(uint)Result_.z]);
            CrossedGlobalEdge[CountCrossedGlobalEdge][1] = make_float2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[((uint)Result_.z + 1) % 3],
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
                                                                 Vertex_Triangle, 1e-7,
                                                                 CrossedGlobalEdge, CountCrossedGlobalEdge, stepNO, i + 1);

            if (Result_.z >= 0)
                goto YK100;
            else
            {

                printf("Warning: I cannot find which edge of element %d does the particle trajectory step over! Particle ID: %d\n", EleID, i + 1);
                goto Debug100;
            }
        };

        /// ----------------crossed edge is identified-----------
        /// ----------------crossed edge is identified-----------
        /// ----------------crossed edge is identified-----------
        uint NumOfElesSharedEdge = EdgesSharedEle_DEV[GlobalEdgeNO].NumSharedEle;

        if (NumOfElesSharedEdge == 1)
        {
            // we have to check if the edge is the target plane or not
            // we have to check if the edge is the target plane or not
            // we have to check if the edge is the target plane or not
            float3 Test3DIntersectionOnCrossedEdge = cuDFNsys::Roate2DPositionTo3D(IntersectionOnCrossedEdge, Frac_DEV[FracID]);
            float *tmp3 = &(Test3DIntersectionOnCrossedEdge.x);
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
            if (stepNO == 1)
            {
                //printf("ParticleID: %d, tmp3[Dir_flow]: %.30f; -1.0f * outletcoordinate: %.30f\n", i + 1, tmp3[Dir_flow], -1.0f * outletcoordinate);
            }
            if (tmp3[Dir_flow] > (-1.0f * outletcoordinate - 1e-4))
            {
                printf("Warning: the particle %d goes above the model!\n", i + 1);
                return;
            }

            // if it is a non-flux plane
            // if it is a non-flux plane
            // if it is a non-flux plane
            if (tmp3[Dir_flow] > outletcoordinate && tmp3[Dir_flow] < (-1.0f * outletcoordinate))
            {
                //printf("Particle %d trajectory meets a non flux bound, before reflection, the trajectory is\n\t%.30f, %.30f\n\t%.30f, %.30f\n", i + 1, InitPos.x, InitPos.y, TargPos.x, TargPos.y);
                TargPos = cuDFNsys::ParticleReflection(TargPos, Vertex_Triangle[(uint)Result_.z], Vertex_Triangle[((uint)Result_.z + 1) % 3]);
                // we do not need to change element ID right now
                // let us determine element ID in next loop
                InitPos = IntersectionOnCrossedEdge;
                //printf("Particle %d trajectory meets a non flux bound, after reflction, the trajectory is\n\t%.30f, %.30f\n\t%.30f, %.30f\n", i + 1, InitPos.x, InitPos.y, TargPos.x, TargPos.y);

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
        }
        else if (NumOfElesSharedEdge == 2)
        {
            /// here, the intersection, if it crosses a vertex, must be very near the vertex and on the bound of the pre element.
            /// here, the intersection, if it crosses a vertex, must be very near the vertex and on the bound of the pre element.
            /// here, the intersection, if it crosses a vertex, must be very near the vertex and on the bound of the pre element.

            bool IfIntersectionPntNearVertex = cuDFNsys::IfParticlePositionNearOneVertexOfElement(IntersectionOnCrossedEdge, Vertex_Triangle, 1e-4);

            if (!IfIntersectionPntNearVertex)
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
                //printf("Particle %d trajectory meets a grid vertex! I won't change element ID. The elementID is %d, the trajectory is\n\t%.30f, %.30f\n\t%.30f, %.30f\n", i + 1, EleID, InitPos.x, InitPos.y, TargPos.x, TargPos.y);
                InitPos = IntersectionOnCrossedEdge;
                continue;
            }
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
            cuDFNsys::WhichElementToGo(EleID,
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
            float3 InitiPos3D = make_float3(IntersectionOnCrossedEdge.x, IntersectionOnCrossedEdge.y, 0.0f);
            float RK_1[3][3];
            Frac_DEV[PreviousFracID].RoationMatrix(RK_1, 23);
            InitiPos3D = cuDFNsys::ProductSquare3Float3(RK_1, InitiPos3D);
            InitiPos3D.x += Frac_DEV[PreviousFracID].Center.x;
            InitiPos3D.y += Frac_DEV[PreviousFracID].Center.y;
            InitiPos3D.z += Frac_DEV[PreviousFracID].Center.z;

            InitiPos3D.x -= Frac_DEV[FracID].Center.x;
            InitiPos3D.y -= Frac_DEV[FracID].Center.y;
            InitiPos3D.z -= Frac_DEV[FracID].Center.z;
            float RK_2[3][3];
            Frac_DEV[FracID].RoationMatrix(RK_2, 32);
            InitiPos3D = cuDFNsys::ProductSquare3Float3(RK_2, InitiPos3D);

            InitPos.x = InitiPos3D.x;
            InitPos.y = InitiPos3D.y;

            // change target position
            // change target position
            // change target position
            float2 tmpDis = make_float2(TargPos.x - IntersectionOnCrossedEdge.x, TargPos.y - IntersectionOnCrossedEdge.y);
            float normDistance = sqrt(tmpDis.x * tmpDis.x + tmpDis.y * tmpDis.y);

            float3 Intersection3D = make_float3(IntersectionOnCrossedEdge.x, IntersectionOnCrossedEdge.y, 0.0f),
                   Target3D = make_float3(TargPos.x, TargPos.y, 0.0f);

            Intersection3D = cuDFNsys::ProductSquare3Float3(RK_1, Intersection3D);
            Intersection3D.x += Frac_DEV[PreviousFracID].Center.x;
            Intersection3D.y += Frac_DEV[PreviousFracID].Center.y;
            Intersection3D.z += Frac_DEV[PreviousFracID].Center.z;
            Target3D = cuDFNsys::ProductSquare3Float3(RK_1, Target3D);
            Target3D.x += Frac_DEV[PreviousFracID].Center.x;
            Target3D.y += Frac_DEV[PreviousFracID].Center.y;
            Target3D.z += Frac_DEV[PreviousFracID].Center.z;

            float3 V_ = make_float3(Target3D.x - Intersection3D.x,
                                    Target3D.y - Intersection3D.y,
                                    Target3D.z - Intersection3D.z);
            float norm_V = sqrt(V_.x * V_.x + V_.y * V_.y + V_.z * V_.z);
            V_.x /= norm_V;
            V_.y /= norm_V;
            V_.z /= norm_V;

            ///// (I -nn) * V
            ///// (I -nn) * V
            ///// (I -nn) * V
            float3 V_new = cuDFNsys::ProjectVToPlaneN(V_, Frac_DEV[FracID].NormalVec);

            //// the angle between the V_new and the 3D velocity vector should be < 90?
            //// get the 3D velocity vector first
            //// the angle between the V_new and the 3D velocity vector should be < 90?
            //// get the 3D velocity vector first
            //// the angle between the V_new and the 3D velocity vector should be < 90?
            //// get the 3D velocity vector first

            float3 NewTarPos3D = make_float3(Intersection3D.x + V_new.x * normDistance - Frac_DEV[FracID].Center.x,
                                             Intersection3D.y + V_new.y * normDistance - Frac_DEV[FracID].Center.y,
                                             Intersection3D.z + V_new.z * normDistance - Frac_DEV[FracID].Center.z);

            NewTarPos3D = cuDFNsys::ProductSquare3Float3(RK_2, NewTarPos3D);
            float2 newTagPos2D = make_float2(NewTarPos3D.x, NewTarPos3D.y);

            // test if the particle goes in an opposite way?
            // test if the particle goes in an opposite way?
            // test if the particle goes in an opposite way?
            float2 EdgeSeg[2];
            EdgeSeg[0].x = Coordinate2D_Vec_dev_ptr[EleID - 1].x[IndexLocal], EdgeSeg[0].y = Coordinate2D_Vec_dev_ptr[EleID - 1].y[IndexLocal];
            EdgeSeg[1].x = Coordinate2D_Vec_dev_ptr[EleID - 1].x[(IndexLocal + 1) % 3], EdgeSeg[1].y = Coordinate2D_Vec_dev_ptr[EleID - 1].y[(IndexLocal + 1) % 3];

            // because FracID changed, the vector of CrossedGlobalEdge should be cleared
            // because FracID changed, the vector of CrossedGlobalEdge should be cleared
            // because FracID changed, the vector of CrossedGlobalEdge should be cleared
            CrossedGlobalEdge[0][0] = EdgeSeg[0];
            CrossedGlobalEdge[0][1] = EdgeSeg[1];

            CountCrossedGlobalEdge = 1;

            float2 CenterGrid;
            CenterGrid.x = (Coordinate2D_Vec_dev_ptr[EleID - 1].x[0] + Coordinate2D_Vec_dev_ptr[EleID - 1].x[1] + Coordinate2D_Vec_dev_ptr[EleID - 1].x[2]) / 3.0f;
            CenterGrid.y = (Coordinate2D_Vec_dev_ptr[EleID - 1].y[0] + Coordinate2D_Vec_dev_ptr[EleID - 1].y[1] + Coordinate2D_Vec_dev_ptr[EleID - 1].y[2]) / 3.0f;

            uint O1 = cuDFNsys::OrientationThree2DPnts(EdgeSeg[0], EdgeSeg[1], newTagPos2D, 1e-7);
            uint O2 = cuDFNsys::OrientationThree2DPnts(EdgeSeg[0], EdgeSeg[1], CenterGrid, 1e-7);

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
            else if (O1 != 0 && O1 == O2)
            {
                TargPos = newTagPos2D;
                whole_Particle_trajectory[0] = InitPos, whole_Particle_trajectory[1] = TargPos;
                //return;

                // float2 Vertex_Triangle_PPP[3];
                // Vertex_Triangle_PPP[0] = make_float2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[0], Coordinate2D_Vec_dev_ptr[EleID - 1].y[0]);
                // Vertex_Triangle_PPP[1] = make_float2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[1], Coordinate2D_Vec_dev_ptr[EleID - 1].y[1]);
                // Vertex_Triangle_PPP[2] = make_float2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[2], Coordinate2D_Vec_dev_ptr[EleID - 1].y[2]);
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
                NewTarPos3D = make_float3(Intersection3D.x - V_new.x * normDistance - Frac_DEV[FracID].Center.x,
                                          Intersection3D.y - V_new.y * normDistance - Frac_DEV[FracID].Center.y,
                                          Intersection3D.z - V_new.z * normDistance - Frac_DEV[FracID].Center.z);
                NewTarPos3D = cuDFNsys::ProductSquare3Float3(RK_2, NewTarPos3D);
                newTagPos2D = make_float2(NewTarPos3D.x, NewTarPos3D.y);
                TargPos = newTagPos2D;
                whole_Particle_trajectory[0] = InitPos, whole_Particle_trajectory[1] = TargPos;
                //return;

                // float2 Vertex_Triangle_PPP[3];
                // Vertex_Triangle_PPP[0] = make_float2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[0], Coordinate2D_Vec_dev_ptr[EleID - 1].y[0]);
                // Vertex_Triangle_PPP[1] = make_float2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[1], Coordinate2D_Vec_dev_ptr[EleID - 1].y[1]);
                // Vertex_Triangle_PPP[2] = make_float2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[2], Coordinate2D_Vec_dev_ptr[EleID - 1].y[2]);
                // printf("O1 != O2, ParticleID %d, eleID: %d;\nPntTrajectory:\n%f, %f\n%f, %f\nTri:\n%f, %f\n%f, %f\n%f, %f\nIntersection3D:\n%f, %f, %f\nV_new:\n%f, %f, %f\nNewTarPos3D:\n%f, %f, %f\n",
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
        };
    };

    if (i + 1 == -1)
    {
    Debug100:;
        {
            printf("Warning: dropped time step!!! ParticleID: %d\nTrajectory:\n\t%.30f, %.30f\n\t%.30f, %.30f\nPre-element ID %d, coordinates:\n\t%.30f, %.30f\n\t%.30f, %.30f\n\t%.30f, %.30f\n", i + 1,
                   designed_particle_trajectory[0].x, designed_particle_trajectory[0].y,
                   designed_particle_trajectory[1].x, designed_particle_trajectory[1].y,
                   InitELeID,
                   Vertex_Triangle_ForVelocity[0].x, Vertex_Triangle_ForVelocity[0].y,
                   Vertex_Triangle_ForVelocity[1].x, Vertex_Triangle_ForVelocity[1].y,
                   Vertex_Triangle_ForVelocity[2].x, Vertex_Triangle_ForVelocity[2].y);
            for (uint k = 0; k < NeighborEleOfOneEle_dev_ptr[InitELeID - 1].NumNeighborEle; ++k)
            {
                uint ele2 = NeighborEleOfOneEle_dev_ptr[InitELeID - 1].EleID[k];

                float2 Vertex_Triangle_PPP[3];
                Vertex_Triangle_PPP[0] = make_float2(Coordinate2D_Vec_dev_ptr[ele2 - 1].x[0], Coordinate2D_Vec_dev_ptr[ele2 - 1].y[0]);
                Vertex_Triangle_PPP[1] = make_float2(Coordinate2D_Vec_dev_ptr[ele2 - 1].x[1], Coordinate2D_Vec_dev_ptr[ele2 - 1].y[1]);
                Vertex_Triangle_PPP[2] = make_float2(Coordinate2D_Vec_dev_ptr[ele2 - 1].x[2], Coordinate2D_Vec_dev_ptr[ele2 - 1].y[2]);
                printf("ParticleID: %d, neighboring elementID: %d:\n\t%.30f, %.30f\n\t%.30f, %.30f\n\t%.30f, %.30f\n",
                       i + 1, ele2,
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
        //     float2 Vertex_Triangle_PPP[3];
        //     Vertex_Triangle_PPP[0] = make_float2(Coordinate2D_Vec_dev_ptr[P_DEV[i].ElementID - 1].x[0], Coordinate2D_Vec_dev_ptr[P_DEV[i].ElementID - 1].y[0]);
        //     Vertex_Triangle_PPP[1] = make_float2(Coordinate2D_Vec_dev_ptr[P_DEV[i].ElementID - 1].x[1], Coordinate2D_Vec_dev_ptr[P_DEV[i].ElementID - 1].y[1]);
        //     Vertex_Triangle_PPP[2] = make_float2(Coordinate2D_Vec_dev_ptr[P_DEV[i].ElementID - 1].x[2], Coordinate2D_Vec_dev_ptr[P_DEV[i].ElementID - 1].y[2]);
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
        //         float2 Vertex_Triangle_PPP[3];
        //         Vertex_Triangle_PPP[0] = make_float2(Coordinate2D_Vec_dev_ptr[ele2 - 1].x[0], Coordinate2D_Vec_dev_ptr[ele2 - 1].y[0]);
        //         Vertex_Triangle_PPP[1] = make_float2(Coordinate2D_Vec_dev_ptr[ele2 - 1].x[1], Coordinate2D_Vec_dev_ptr[ele2 - 1].y[1]);
        //         Vertex_Triangle_PPP[2] = make_float2(Coordinate2D_Vec_dev_ptr[ele2 - 1].x[2], Coordinate2D_Vec_dev_ptr[ele2 - 1].y[2]);
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
}; // namespace cuDFNsys
