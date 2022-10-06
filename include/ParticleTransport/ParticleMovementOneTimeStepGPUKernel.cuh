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
#include "../Mesh/EleCoor.cuh"
#include "../RandomFunction/RandomFunction.cuh"
#include "AngleBetweenTwoNeighboringTriangles.cuh"
#include "EdgeToEle.cuh"
#include "IdentifyParticleCrossesWhichEdge.cuh"
#include "IfParticleOnBoundOfElement.cuh"
#include "IfParticlePositionInNeighboringElement.cuh"
#include "IfParticlePositionNearOneVertexOfElement.cuh"
#include "Particle.cuh"
#include "ParticleReflection.cuh"
#include "Roate2DPositionTo3D.cuh"
#include "RotationRemainningTrajectoryBetweenTwo3DTriangles.cuh"
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
                                                     //cuDFNsys::NeighborEle *NeighborEleOfOneEle_dev_ptr,
                                                     uint *EleToFracID_ptr,
                                                     T *velocity_ptr,
                                                     uint Dir_flow,
                                                     T outletcoordinate,
                                                     int count,
                                                     int numElements,
                                                     uint stepNO,
                                                     uint *Particle_runtime_error_dev_pnt)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= count)
        return;

    /// ----------------------- debug -----------------------
    /// ----------------------- debug -----------------------
    /// ----------------------- debug -----------------------
    // zP_DEV[i].ElementID = 11847;
    // P_DEV[i].Position2D = cuDFNsys::MakeVector2<T>(33.9676100000000005252331902738660573959351, -5.3694499999999996120436662749852985143661);
    /// ----------------------- debug -----------------------
    /// ----------------------- debug -----------------------
    /// ----------------------- debug -----------------------

    if (P_DEV[i].ParticleID == -1)
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

    /// the dimension of Veloc_p is 2, and the unit is L^2 * S^-1
    /// we have to convert it to L * S^-1
    T conductivity_k = Frac_DEV[FracID].Conductivity;
    T b_aperture = pow(conductivity_k * 12.0, 1.0 / 3.0);
    Veloc_p.x /= b_aperture, Veloc_p.y /= b_aperture;

    /// ----------------------- debug -----------------------
    /// ----------------------- debug -----------------------
    /// ----------------------- debug -----------------------
    /// bool IfTargPosInGrid_yy = cuDFNsys::IfPntInside2DConvexPoly<T>(InitPos, Vertex_Triangle_ForVelocity, 3);
    /// if (!IfTargPosInGrid_yy)
    /// {
    ///     printf("not inside~\n");
    ///     return;
    /// }
    //
    //printf("velocity 2D: %.40f, %.40f\n", Veloc_p.x, Veloc_p.y);

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
    /// TargPos = cuDFNsys::MakeVector2<T>(34.0291899999999998271960066631436347961426, -5.3011900000000000687805368215776979923248);
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

    cuDFNsys::Vector2<T> CrossedGlobalEdge[_SizeOfArray_CrossedGlobalEdge_][2];
    int CountCrossedGlobalEdge = 0;

    cuDFNsys::Vector2<T> whole_Particle_trajectory[2] = {InitPos, TargPos};

    for (uint Loop_time = 1;; Loop_time++)
    {
        /// printf("\n%%--------------------\n%% The looptime: %d, fracid: %d, elementid: %d\n %%Present trajectory:\n%.40f, %.40f\n%.40f, %.40f,\n\n",
        ///        Loop_time, FracID, EleID, InitPos.x, InitPos.y, TargPos.x, TargPos.y);
        /// ////// debug: turn trajectory to 3D --------------------
        /// ////// debug: turn trajectory to 3D --------------------
        /// ////// debug: turn trajectory to 3D --------------------
        /// if (1)
        /// {
        ///     cuDFNsys::Vector3<T> InitiPos3D = cuDFNsys::MakeVector3(InitPos.x, InitPos.y, (T)0.0);
        ///     cuDFNsys::Vector3<T> TargPos3D = cuDFNsys::MakeVector3(TargPos.x, TargPos.y, (T)0.0);
        ///     T RK_1[3][3];
        ///     Frac_DEV[FracID].RoationMatrix(RK_1, 23);
        ///     InitiPos3D = cuDFNsys::ProductSquare3Vector3<T>(RK_1, InitiPos3D);
        ///     TargPos3D = cuDFNsys::ProductSquare3Vector3<T>(RK_1, TargPos3D);
        ///     InitiPos3D.x += Frac_DEV[FracID].Center.x;
        ///     InitiPos3D.y += Frac_DEV[FracID].Center.y;
        ///     InitiPos3D.z += Frac_DEV[FracID].Center.z;
        ///     TargPos3D.x += Frac_DEV[FracID].Center.x;
        ///     TargPos3D.y += Frac_DEV[FracID].Center.y;
        ///     TargPos3D.z += Frac_DEV[FracID].Center.z;
        ///     printf("Trajectory3D:\n\t%.40f, %.40f, %.40f\n\t%.40f, %.40f, %.40f\n",
        ///            InitiPos3D.x,
        ///            InitiPos3D.y,
        ///            InitiPos3D.z,
        ///            TargPos3D.x,
        ///            TargPos3D.y,
        ///            TargPos3D.z);
        /// }

        if (Loop_time >= _ParTran_MaxLoopTimes || CountCrossedGlobalEdge == _SizeOfArray_CrossedGlobalEdge_ - 1)
        {
            printf("Particle %d, Loop times is too large!\n", i + 1);
            goto Debug100;
        }
        cuDFNsys::Vector2<T> Vertex_Triangle[3];
        Vertex_Triangle[0] = cuDFNsys::MakeVector2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[0], Coordinate2D_Vec_dev_ptr[EleID - 1].y[0]);
        Vertex_Triangle[1] = cuDFNsys::MakeVector2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[1], Coordinate2D_Vec_dev_ptr[EleID - 1].y[1]);
        Vertex_Triangle[2] = cuDFNsys::MakeVector2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[2], Coordinate2D_Vec_dev_ptr[EleID - 1].y[2]);

        // printf("Element id: %d, 2D coordinates:\n\t%.40f, %.40f\n\t%.40f, %.40f\n\t%.40f, %.40f\n",
        //        EleID, Vertex_Triangle[0].x, Vertex_Triangle[0].y, Vertex_Triangle[1].x, Vertex_Triangle[1].y, Vertex_Triangle[2].x, Vertex_Triangle[2].y);
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
                                                                         (T)1e-12);

        //printf("IfTargPosOnGridBound: %d\n", IfTargPosOnGridBound);
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
                P_DEV[i].ParticleID = -1; // reaches outlet plane
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
        cuDFNsys::Vector3<T> Result_ = cuDFNsys::IdentifyParticleCrossesWhichEdge(P_trajectory_, Vertex_Triangle, (T)1e-12, CrossedGlobalEdge, CountCrossedGlobalEdge, stepNO, i + 1);

        uint GlobalEdgeNO = 0;
        cuDFNsys::Vector2<T> IntersectionOnCrossedEdge;
        if (Result_.z >= 0)
        {
        YK100:;
            IntersectionOnCrossedEdge = cuDFNsys::MakeVector2(Result_.x, Result_.y);
            uint *tmp_k = &(EdgeNO.x);
            GlobalEdgeNO = tmp_k[(uint)Result_.z];
            // printf("\nResult_.z >= 0; Result_.z = %.40f, edgeno: %d, %d, %d\nIntersection2D: %.40f, %.40f\n",
            //        Result_.z, EdgeNO.x, EdgeNO.y, EdgeNO.z,
            //        IntersectionOnCrossedEdge.x, IntersectionOnCrossedEdge.y);

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
                /// in case of: a particle went through a fracture trace and got into anther fracture, and the rotation
                /// made the initial point out of the element, because the previous trajectory crossed a vertex!!!

                bool IfInitPosInGrid_yy = cuDFNsys::IfPntInside2DConvexPoly<T>(InitPos, Vertex_Triangle, 3);
                bool IfTargPosInGrid_yy = cuDFNsys::IfPntInside2DConvexPoly<T>(TargPos, Vertex_Triangle, 3);

                //printf("IfInitPosInGrid_yy %d, IfTargPosInGrid_yy %d\n", IfInitPosInGrid_yy, IfTargPosInGrid_yy);

                if (!IfInitPosInGrid_yy && !IfTargPosInGrid_yy)
                    for (uint yi = 0; yi < 3; ++yi)
                    {
                        cuDFNsys::Vector2<T> p2, q2;
                        p2 = Vertex_Triangle[yi];
                        q2 = Vertex_Triangle[(yi + 1) % 3];

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

                        /// if p2 or q2 is very close to InitPos ?
                        bool close_p2_or_q2 = false;

                        T dis1 = sqrt((p2.x - InitPos.x) * (p2.x - InitPos.x) + (p2.y - InitPos.y) * (p2.y - InitPos.y));
                        T dis2 = sqrt((q2.x - InitPos.x) * (q2.x - InitPos.x) + (q2.y - InitPos.y) * (q2.y - InitPos.y));

                        if (dis1 < 1e-7 || dis2 < 1e-7)
                            close_p2_or_q2 = true;

                        // printf("dis1 %.40f, dis2 %.40f\n", dis1, dis2);
                        if (close_p2_or_q2)
                        {
                            IntersectionOnCrossedEdge = InitPos;
                            uint *tmp_k = &(EdgeNO.x);
                            GlobalEdgeNO = tmp_k[yi];
                            // printf("\nResult_.z >= 0; Result_.z = %.40f, edgeno: %d, %d, %d\nIntersection2D: %.40f, %.40f\n",
                            //        Result_.z, EdgeNO.x, EdgeNO.y, EdgeNO.z,
                            //        IntersectionOnCrossedEdge.x, IntersectionOnCrossedEdge.y);

                            CrossedGlobalEdge[CountCrossedGlobalEdge][0] = cuDFNsys::MakeVector2(p2.x,
                                                                                                 p2.y);
                            CrossedGlobalEdge[CountCrossedGlobalEdge][1] = cuDFNsys::MakeVector2(q2.x,
                                                                                                 q2.y);
                            CountCrossedGlobalEdge++;

                            break;
                        }
                        else if (!close_p2_or_q2 && yi == 2)
                        {
                            /////
                            printf("Warning: I cannot find which edge of element %d does the particle trajectory step over I! Particle ID: %d\n", EleID, i + 1);
                            goto Debug100;
                        }
                    }
                else
                {
                    printf("Warning: I cannot find which edge of element %d does the particle trajectory step over II! Particle ID: %d\n", EleID, i + 1);
                    goto Debug100;
                }
            }
        };

        /// ---------------- if the intersection point is on the inlet plane, and the target point is above the inlet plane? ---------
        /// ---------------- if the intersection point is on the inlet plane, and the target point is above the inlet plane? --------
        /// ---------------- if the intersection point is on the inlet plane, and the target point is above the inlet plane? --------
        /// ---------------- this happens, particularly for the first step when the molecular diffusion is large--------------
        /// ---------------- this happens, particularly for the first step when the molecular diffusion is large--------------
        /// ---------------- this happens, particularly for the first step when the molecular diffusion is large--------------
        cuDFNsys::Vector3<T> IntersectionPnt_3D = cuDFNsys::Roate2DPositionTo3D<T>(IntersectionOnCrossedEdge,
                                                                                   Frac_DEV[EleToFracID_ptr[EleID - 1]]);

        cuDFNsys::Vector3<T> TargetPnt_3D = cuDFNsys::Roate2DPositionTo3D<T>(TargPos,
                                                                             Frac_DEV[EleToFracID_ptr[EleID - 1]]);
        T *tmp4 = &(IntersectionPnt_3D.x);
        if (abs(tmp4[Dir_flow] - ((T)-1.0) * outletcoordinate) < 1e-4)
        {
            tmp4 = &(TargetPnt_3D.x);

            //printf("abs(tmp4[Dir_flow] - ((T)-1.0) * outletcoordinate): %.40f, %.40f, %.40f, %d\n", abs(tmp4[Dir_flow] - ((T)-1.0) * outletcoordinate), tmp4[Dir_flow], ((T)-1.0) * outletcoordinate, EleToFracID_ptr[P_DEV[i].ElementID - 1]);

            if (tmp4[Dir_flow] < ((T)-1.0) * outletcoordinate && tmp4[Dir_flow] > outletcoordinate)
            {
                // do nothing
            }
            else
            {
                printf("Warning: the particle %d goes above the model or moves along the inlet plane!\n", i + 1);
                // delete this particle in the future
                // goto Debug100;
                return;
            }
        }

        /// ----------------------- if the intersection point is on the outlet plane
        /// ----------------------- if the intersection point is on the outlet plane
        /// ----------------------- if the intersection point is on the outlet plane
        tmp4 = &(IntersectionPnt_3D.x);
        if (abs(tmp4[Dir_flow] - outletcoordinate) < 1e-4)
        {
            P_DEV[i].ParticleID = -1;
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
                                                                                   Frac2Dvertex, Frac_DEV[FracID].NumVertsTruncated,
                                                                                   (T)1e-12,
                                                                                   &edgeNO);
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
            //printf("~~\n");
            continue;
        }

        delete[] Frac2Dvertex;
        Frac2Dvertex = NULL;

        /// ----------------crossed edge is identified-----------
        /// ----------------crossed edge is identified-----------
        /// ----------------crossed edge is identified-----------
        //printf("\n\nNumOfElesSharedEdge = %d\n\n", NumOfElesSharedEdge);

        if (NumOfElesSharedEdge == 1)
        {
            // we have to check if the edge is the target plane or not
            // we have to check if the edge is the target plane or not
            // we have to check if the edge is the target plane or not
            cuDFNsys::Vector3<T> Test3DIntersectionOnCrossedEdge = cuDFNsys::Roate2DPositionTo3D<T>(IntersectionOnCrossedEdge, Frac_DEV[FracID]);
            T *tmp3 = &(Test3DIntersectionOnCrossedEdge.x);
            if (abs(tmp3[Dir_flow] - outletcoordinate) < 1e-4)
            {
                P_DEV[i].ParticleID = -1;
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

                    // printf("Particle %d trajectory meets a non flux bound, before reflection, the trajectory is\n\t%.40f, %.40f\n\t%.40f, %.40f\n",
                    //        i + 1, InitPos.x, InitPos.y, TargPos.x, TargPos.y);
                    TargPos = cuDFNsys::ParticleReflection<T>(TargPos, Vertex_Triangle[(uint)Result_.z], Vertex_Triangle[((uint)Result_.z + 1) % 3]);
                    // we do not need to change element ID right now
                    // let us determine element ID in next loop
                    InitPos = IntersectionOnCrossedEdge;
                    // printf("Particle %d trajectory meets a non flux bound, after reflction, the trajectory is\n\t%.40f, %.40f\n\t%.40f, %.40f\n", i + 1, InitPos.x, InitPos.y, TargPos.x, TargPos.y);

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

            //printf("Particle %d trajectory meets a grid edge, the elementID is %d, the cross global edgeNO is %d, the trajectory is\n\t%.40f, %.40f\n\t%.40f, %.40f\n", i + 1, EleID, EleID * 3 - 3 + (uint)Result_.z, InitPos.x, InitPos.y, TargPos.x, TargPos.y);
            EleID = (EdgesSharedEle_DEV[GlobalEdgeNO].EleID[0] == EleID ? EdgesSharedEle_DEV[GlobalEdgeNO].EleID[1] : EdgesSharedEle_DEV[GlobalEdgeNO].EleID[0]);
            //P_DEV[i].ElementID = EleID;
            InitPos = IntersectionOnCrossedEdge;
            //printf("Particle %d trajectory meets a grid edge, the elementID change to %d, the trajectory is\n\t%.40f, %.40f\n\t%.40f, %.40f\n\tAccumulated crossed edge global NO:\n", i + 1, EleID, InitPos.x, InitPos.y, TargPos.x, TargPos.y);
            FracID = EleToFracID_ptr[EleID - 1];                              // from 0
            EdgeNO = make_uint3(EleID * 3 - 3, EleID * 3 - 2, EleID * 3 - 1); // from 0

            //for (uint y = 0; y < CountCrossedGlobalEdge; ++y)
            //printf("\tParticle %d, crossed global edgeNO:\n\t\t%.40f, %.40f\n\t\t%.40f, %.40f\n", i + 1, CrossedGlobalEdge[y][0].x, CrossedGlobalEdge[y][0].y, CrossedGlobalEdge[y][1].x, CrossedGlobalEdge[y][1].y);
            continue;
        }
        else
        {
            // may go to another fracture
            // may go to another fracture
            // may go to another fracture

            uint PreviousFracID = FracID;
            uint PreEleID = EleID;
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
                                          ifAllsharedEdgeVelocityPositive);
            // printf("ifAllsharedEdgeVelocityPositive: %d\n", ifAllsharedEdgeVelocityPositive);

            if (ifAllsharedEdgeVelocityPositive == false && (newELEID_ == -1 || newFracID_ == -1 || IndexLocal == -1))
            {
                printf("Warning: illeagal indice appear, because the determination of which element to go is failed!\nnewELEID_: %d\nnewFracID_: %d\nIndexLocal: %d\n",
                       newELEID_,
                       newFracID_,
                       IndexLocal);
                goto Debug100;
            }

            if (ifAllsharedEdgeVelocityPositive == true)
            {
                /// the particle trajectory goes back to the previous element, by reflection
                /// the particle trajectory goes back to the previous element, by reflection
                /// the particle trajectory goes back to the previous element, by reflection

                cuDFNsys::Vector2<T> EdgeSeg_II[2];
                EdgeSeg_II[0].x = Coordinate2D_Vec_dev_ptr[EleID - 1].x[(uint)Result_.z];
                EdgeSeg_II[0].y = Coordinate2D_Vec_dev_ptr[EleID - 1].y[(uint)Result_.z];
                EdgeSeg_II[1].x = Coordinate2D_Vec_dev_ptr[EleID - 1].x[((uint)Result_.z + 1) % 3];
                EdgeSeg_II[1].y = Coordinate2D_Vec_dev_ptr[EleID - 1].y[((uint)Result_.z + 1) % 3];

                InitPos = IntersectionOnCrossedEdge;
                TargPos = cuDFNsys::ParticleReflection<T>(TargPos, EdgeSeg_II[0], EdgeSeg_II[1]);

                CrossedGlobalEdge[0][0] = EdgeSeg_II[0];
                CrossedGlobalEdge[0][1] = EdgeSeg_II[1];

                CountCrossedGlobalEdge = 1;

                continue;
            };

            FracID = (uint)newFracID_;
            EleID = (uint)newELEID_;
            EdgeNO = make_uint3(EleID * 3 - 3, EleID * 3 - 2, EleID * 3 - 1); // from 0

            if (FracID == PreviousFracID) // did not change fracture ID
            {
                InitPos = IntersectionOnCrossedEdge;
                continue;
            }

            //-------------
            //-------------
            //-------------
            uint designed_next_EleID = 0;
            // if there is no intersection of fractures, then the particle must get into this element!
            uint localedgeno_ppp = 0;

            //printf("PreEleID: %d, PreviousFracID: %d\n", PreEleID, PreviousFracID);
            for (uint p = 0; p < NumOfElesSharedEdge; ++p)
            {
                //printf("EdgesSharedEle_DEV[GlobalEdgeNO].EleID[p]: %d, fracid: %d\n", EdgesSharedEle_DEV[GlobalEdgeNO].EleID[p], EleToFracID_ptr[EdgesSharedEle_DEV[GlobalEdgeNO].EleID[p] - 1]);
                if (EdgesSharedEle_DEV[GlobalEdgeNO].EleID[p] != PreEleID && EleToFracID_ptr[EdgesSharedEle_DEV[GlobalEdgeNO].EleID[p] - 1] == PreviousFracID)
                {
                    designed_next_EleID = EdgesSharedEle_DEV[GlobalEdgeNO].EleID[p];
                    localedgeno_ppp = EdgesSharedEle_DEV[GlobalEdgeNO].LocalEdgeNO[p];
                    break;
                }
                else if (p == NumOfElesSharedEdge - 1)
                {
                    printf("warning: I did not find the neighboring element which is on the same plane of previous element!\n");
                    goto Debug100;
                };
            }

            T RK_1[3][3];
            Frac_DEV[PreviousFracID].RoationMatrix(RK_1, 23);

            T RK_2[3][3];
            Frac_DEV[FracID].RoationMatrix(RK_2, 32);

            T RK_3[3][3];
            Frac_DEV[FracID].RoationMatrix(RK_3, 23);

            cuDFNsys::Vector3<T> designed_element[3], next_element[3];
            for (uint yr = 0; yr < 3; ++yr)
            {
                cuDFNsys::Vector3<T> Pnt_tt = cuDFNsys::MakeVector3(Coordinate2D_Vec_dev_ptr[designed_next_EleID - 1].x[yr],
                                                                    Coordinate2D_Vec_dev_ptr[designed_next_EleID - 1].y[yr],
                                                                    (T)0.);
                Pnt_tt = cuDFNsys::ProductSquare3Vector3<T>(RK_1, Pnt_tt);
                designed_element[yr] = cuDFNsys::MakeVector3(Pnt_tt.x + Frac_DEV[PreviousFracID].Center.x,
                                                             Pnt_tt.y + Frac_DEV[PreviousFracID].Center.y,
                                                             Pnt_tt.z + Frac_DEV[PreviousFracID].Center.z);
            }

            for (uint yr = 0; yr < 3; ++yr)
            {
                cuDFNsys::Vector3<T> Pnt_tt = cuDFNsys::MakeVector3(Coordinate2D_Vec_dev_ptr[EleID - 1].x[yr],
                                                                    Coordinate2D_Vec_dev_ptr[EleID - 1].y[yr],
                                                                    (T)0.);
                Pnt_tt = cuDFNsys::ProductSquare3Vector3<T>(RK_3, Pnt_tt);
                next_element[yr] = cuDFNsys::MakeVector3(Pnt_tt.x + Frac_DEV[FracID].Center.x,
                                                         Pnt_tt.y + Frac_DEV[FracID].Center.y,
                                                         Pnt_tt.z + Frac_DEV[FracID].Center.z);
            }

            ///------------------ d1 is a directional vector on the true next element and perpendicular to rotation axis
            ///------------------ d2 is a directional vector on the designed next element and perpendicular to rotation axis
            // cuDFNsys::Vector3<T> d1, d2;
            //----------------------------------------------------------------------------------------------------------
            // printf("PreEleID: %d, EleID: %d, designed_next_EleID: %d\n", PreEleID, EleID, designed_next_EleID);
            T angle_ = cuDFNsys::AngleBetweenTwoNeighboringTriangles<T>(designed_element, next_element,
                                                                        localedgeno_ppp, IndexLocal /*, d1, d2*/);

            cuDFNsys::Vector3<T> Target3D = cuDFNsys::MakeVector3(TargPos.x, TargPos.y, (T)0.0f);
            Target3D = cuDFNsys::ProductSquare3Vector3<T>(RK_1, Target3D);
            Target3D.x += Frac_DEV[PreviousFracID].Center.x;
            Target3D.y += Frac_DEV[PreviousFracID].Center.y;
            Target3D.z += Frac_DEV[PreviousFracID].Center.z; /// rotate this point

            cuDFNsys::Vector3<T> designed_element_ppp[3];

            designed_element_ppp[0] = Target3D;
            designed_element_ppp[1] = designed_element[localedgeno_ppp];
            designed_element_ppp[2] = designed_element[(localedgeno_ppp + 1) % 3];

            cuDFNsys::Vector3<T> new_Target3D;
            //printf("angle_: %.40f\n", angle_);
            T angle_ppp[2];
            for (uint yr = 0; yr < 2; ++yr)
            {
                new_Target3D = cuDFNsys::RotationRemainningTrajectoryBetweenTwo3DTriangles<T>(designed_element_ppp, angle_);
                // printf("new_Target3D: %.40f, %.40f, %.40f\n", new_Target3D.x, new_Target3D.y, new_Target3D.z);
                //------------test--------
                //cuDFNsys::Vector3<T> d3, d4;

                /// test if the new_Target3D is in the correct 3D element?
                /// test if the new_Target3D is in the correct 3D element?
                /// test if the new_Target3D is in the correct 3D element?

                cuDFNsys::Vector3<T> next_element_ppp[3];
                next_element_ppp[0] = new_Target3D;
                next_element_ppp[1] = next_element[IndexLocal];
                next_element_ppp[2] = next_element[(IndexLocal + 1) % 3];

                angle_ppp[yr] = cuDFNsys::AngleBetweenTwoNeighboringTriangles<T>(next_element_ppp, next_element, 1,
                                                                                 IndexLocal /*, d3, d4*/);

                //printf("angle_ppp[%d]: %.40f\n", yr, angle_ppp[yr]);
                T ang_degree = angle_ppp[yr] * 180.0 / M_PI;

                T error_td = 0.3;

                if (ang_degree > error_td && yr < 1)
                {
                    angle_ = (T)2.0 * M_PI - angle_;
                    continue;
                }
                else if (ang_degree > error_td && yr == 1)
                {
                    printf("warning: I cannot correctly rotate the remainning trajectory when particle goes through a fracture trace! angle_ppp: %.40f, %.40f\n",
                           angle_ppp[0] * 180.0f / M_PI, angle_ppp[1] * 180.0f / M_PI);
                    goto Debug100;
                }
                else
                    break;
            };

            new_Target3D = cuDFNsys::MakeVector3(new_Target3D.x - Frac_DEV[FracID].Center.x,
                                                 new_Target3D.y - Frac_DEV[FracID].Center.y,
                                                 new_Target3D.z - Frac_DEV[FracID].Center.z);

            new_Target3D = cuDFNsys::ProductSquare3Vector3<T>(RK_2, new_Target3D);

            TargPos = cuDFNsys::MakeVector2(new_Target3D.x, new_Target3D.y);

            cuDFNsys::Vector3<T> InitiPos3D = cuDFNsys::MakeVector3(IntersectionOnCrossedEdge.x, IntersectionOnCrossedEdge.y, (T)0.0f);
            InitiPos3D = cuDFNsys::ProductSquare3Vector3<T>(RK_1, InitiPos3D);
            InitiPos3D.x += Frac_DEV[PreviousFracID].Center.x;
            InitiPos3D.y += Frac_DEV[PreviousFracID].Center.y;
            InitiPos3D.z += Frac_DEV[PreviousFracID].Center.z;
            InitiPos3D.x -= Frac_DEV[FracID].Center.x;
            InitiPos3D.y -= Frac_DEV[FracID].Center.y;
            InitiPos3D.z -= Frac_DEV[FracID].Center.z;
            InitiPos3D = cuDFNsys::ProductSquare3Vector3<T>(RK_2, InitiPos3D);
            InitPos.x = InitiPos3D.x;
            InitPos.y = InitiPos3D.y;
            whole_Particle_trajectory[0] = InitPos, whole_Particle_trajectory[1] = TargPos;

            // because FracID changed, the vector of CrossedGlobalEdge should be cleared
            // because FracID changed, the vector of CrossedGlobalEdge should be cleared
            // because FracID changed, the vector of CrossedGlobalEdge should be cleared
            cuDFNsys::Vector2<T> EdgeSeg[2];
            EdgeSeg[0].x = Coordinate2D_Vec_dev_ptr[EleID - 1].x[IndexLocal],
            EdgeSeg[0].y = Coordinate2D_Vec_dev_ptr[EleID - 1].y[IndexLocal];
            EdgeSeg[1].x = Coordinate2D_Vec_dev_ptr[EleID - 1].x[(IndexLocal + 1) % 3],
            EdgeSeg[1].y = Coordinate2D_Vec_dev_ptr[EleID - 1].y[(IndexLocal + 1) % 3];

            CrossedGlobalEdge[0][0] = EdgeSeg[0];
            CrossedGlobalEdge[0][1] = EdgeSeg[1];

            CountCrossedGlobalEdge = 1;

            continue;
        };
    };

    if (i + 1 == -1)
    {
    Debug100:;
        {
            printf("Warning: dropped time step!!! ParticleID: %d\n%%Trajectory:\n\t%.40f, %.40f\n\t%.40f, %.40f\nPre-element ID %d, fracID: %d, coordinates:\n\t%.40f, %.40f\n\t%.40f, %.40f\n\t%.40f, %.40f\n", i + 1,
                   designed_particle_trajectory[0].x, designed_particle_trajectory[0].y,
                   designed_particle_trajectory[1].x, designed_particle_trajectory[1].y,
                   InitELeID, EleToFracID_ptr[InitELeID - 1],
                   Vertex_Triangle_ForVelocity[0].x, Vertex_Triangle_ForVelocity[0].y,
                   Vertex_Triangle_ForVelocity[1].x, Vertex_Triangle_ForVelocity[1].y,
                   Vertex_Triangle_ForVelocity[2].x, Vertex_Triangle_ForVelocity[2].y);

            //for (uint k = 0; k < NeighborEleOfOneEle_dev_ptr[InitELeID - 1].NumNeighborEle; ++k)
            //{
            //    uint ele2 = NeighborEleOfOneEle_dev_ptr[InitELeID - 1].EleID[k];
            //    cuDFNsys::Vector2<T> Vertex_Triangle_PPP[3];
            //    Vertex_Triangle_PPP[0] = cuDFNsys::MakeVector2(Coordinate2D_Vec_dev_ptr[ele2 - 1].x[0], Coordinate2D_Vec_dev_ptr[ele2 - 1].y[0]);
            //    Vertex_Triangle_PPP[1] = cuDFNsys::MakeVector2(Coordinate2D_Vec_dev_ptr[ele2 - 1].x[1], Coordinate2D_Vec_dev_ptr[ele2 - 1].y[1]);
            //    Vertex_Triangle_PPP[2] = cuDFNsys::MakeVector2(Coordinate2D_Vec_dev_ptr[ele2 - 1].x[2], Coordinate2D_Vec_dev_ptr[ele2 - 1].y[2]);
            //    printf("%% ParticleID: %d, neighboring elementID: %d, fracID: %d:\n\t%.40f, %.40f\n\t%.40f, %.40f\n\t%.40f, %.40f\n",
            //           i + 1, ele2, EleToFracID_ptr[ele2 - 1],
            //           Vertex_Triangle_PPP[0].x, Vertex_Triangle_PPP[0].y,
            //           Vertex_Triangle_PPP[1].x, Vertex_Triangle_PPP[1].y,
            //           Vertex_Triangle_PPP[2].x, Vertex_Triangle_PPP[2].y);
            //}
        };
        Particle_runtime_error_dev_pnt[i] = 1;
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
        //     printf("Time step infomation. ParticleID: %d\nTrajectory:\n\t%.40f, %.40f\n\t%.40f, %.40f\nPre-element ID %d, coordinates:\n\t%.40f, %.40f\n\t%.40f, %.40f\n\t%.40f, %.40f\nParticleID: %d, Now, element ID: %d, coordinates:\n\t%.40f, %.40f\n\t%.40f, %.40f\n\t%.40f, %.40f\n",
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
        //         printf("ParticleID: %d, neighboring elementID: %d:\n\t%.40f, %.40f\n\t%.40f, %.40f\n\t%.40f, %.40f\n",
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
                                                                      //cuDFNsys::NeighborEle *NeighborEleOfOneEle_dev_ptr,
                                                                      uint *EleToFracID_ptr,
                                                                      double *velocity_ptr,
                                                                      uint Dir_flow,
                                                                      double outletcoordinate,
                                                                      int count,
                                                                      int numElements,
                                                                      uint stepNO, uint *Particle_runtime_error_dev_pnt);
template __global__ void ParticleMovementOneTimeStepGPUKernel<float>(unsigned long seed,
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
                                                                     uint stepNO, uint *Particle_runtime_error_dev_pnt);
}; // namespace cuDFNsys
