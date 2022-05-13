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
#include "IdentifyParticleElementAllElements.cuh"
#include "IfParticleInFrac.cuh"
#include "IfParticleInsideGrid.cuh"
#include "Particle.cuh"
#include "ParticleReflection.cuh"
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
                                                     int numElements)
{
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

    // ----------------check if the new position is in the grid----------------
    // ----------------check if the new position is in the grid----------------
    // ----------------check if the new position is in the grid----------------

    uint Record_1stEleID = EleID;
    uint count_loop = 0;

StartOfCheckParticlePosition:;

    count_loop++;

    float2 Vertex_Triangle[3];
    Vertex_Triangle[0] = make_float2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[0], Coordinate2D_Vec_dev_ptr[EleID - 1].y[0]);
    Vertex_Triangle[1] = make_float2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[1], Coordinate2D_Vec_dev_ptr[EleID - 1].y[1]);
    Vertex_Triangle[2] = make_float2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[2], Coordinate2D_Vec_dev_ptr[EleID - 1].y[2]);

    // bool IfInGrid_3 = cuDFNsys::IfParticleInsideGrid(InitPos,
    //                                                  Vertex_Triangle);
    // if (IfInGrid_3 == false)
    // {
    //     printf("IfInGrid_3: false, particle %d, eleID: %d\nPnt:\n%f, %f\nTri:\n%f, %f\n%f, %f\n%f, %f\n\n", i + 1, EleID,
    //            InitPos.x, InitPos.y,
    //            Vertex_Triangle[0].x, Vertex_Triangle[0].y,
    //            Vertex_Triangle[1].x, Vertex_Triangle[1].y,
    //            Vertex_Triangle[2].x, Vertex_Triangle[2].y);
    // }

    bool IfInGrid = cuDFNsys::IfParticleInsideGrid(TargPos,
                                                   Vertex_Triangle);

    //-----------// in the grid
    //-----------// in the grid
    //-----------// in the grid
    if (IfInGrid == true)
    {
        P_DEV[i].Position2D = TargPos;
        return;
    }
    //-------- out of the grid---------
    //-------- out of the grid---------
    //-------- out of the grid---------
    float2 ParticleTrajectory[2] = {InitPos, TargPos};

    uint init_checking_edgeNO = (uint)floor(cuDFNsys::RandomUniform(0.0f, 2.99f, curand_uniform(&state)));

    float3 ResultParticle = cuDFNsys::IdentifyParticleCrossesWhichEdge(ParticleTrajectory, Vertex_Triangle, init_checking_edgeNO, 1e-5);
    uint CrossedGlobalEdgeNO = 0;

    float2 IntersectionPnt = make_float2(ResultParticle.x, ResultParticle.y);

    if (count_loop >= 10)
    {
        printf("Warning: particle trajectory is too large?\nParticleID %d, eleID: %d; Record_1stEleID: %d\nPntTrajectory:\n%f, %f\n%f, %f\nTri:\n%f, %f\n%f, %f\n%f, %f\nIntersection:\n%f, %f\nEdgeNO: %f\n\n", i + 1, EleID, Record_1stEleID,
               InitPos.x, InitPos.y, TargPos.x, TargPos.y,
               Vertex_Triangle[0].x, Vertex_Triangle[0].y,
               Vertex_Triangle[1].x, Vertex_Triangle[1].y,
               Vertex_Triangle[2].x, Vertex_Triangle[2].y,
               IntersectionPnt.x, IntersectionPnt.y, ResultParticle.z);
    };

    if (ResultParticle.z >= 0)
    {
        uint *tmp_k = &(EdgeNO.x);
        CrossedGlobalEdgeNO = tmp_k[(uint)ResultParticle.z];

        // printf("P_ID %d, CrossedGlobalEdgeNO: %d, intersection:\n%f, %f\nsegment:\n%f, %f\n%f, %f\nTri:%f, %f\n%f, %f\n%f, %f\n\n", i + 1, CrossedGlobalEdgeNO,
        //        IntersectionPnt.x, IntersectionPnt.y,
        //        InitPos.x, InitPos.y, TargPos.x, TargPos.y,
        //        Vertex_Triangle[0].x, Vertex_Triangle[0].y,
        //        Vertex_Triangle[1].x, Vertex_Triangle[1].y,
        //        Vertex_Triangle[2].x, Vertex_Triangle[2].y);
    }
    else
    {
        // bool IfInGrid_2 = cuDFNsys::IfParticleInsideGrid(InitPos,
        //                                                  Vertex_Triangle);
        // printf("Warning: I cannot identify which edge of the grid does the particle %d cross through\nSo I just drop this time step!\nThe triangle:\n%f, %f\n%f, %f\n%f, %f\nThe segment:\n%f, %f\n%f, %f\nIfInGrid_2: %d\n\n",
        //        i + 1,
        //        Vertex_Triangle[0].x, Vertex_Triangle[0].y,
        //        Vertex_Triangle[1].x, Vertex_Triangle[1].y,
        //        Vertex_Triangle[2].x, Vertex_Triangle[2].y,
        //        InitPos.x, InitPos.y,
        //        TargPos.x, TargPos.y, IfInGrid_2);
        // return;

        bool GGH = cuDFNsys::IdentifyParticleElementAllElements(InitPos, TargPos, EleID, numElements, EleToFracID_ptr,
                                                                FracID, NeighborEleOfOneEle_dev_ptr[P_DEV[i].ElementID - 1], Coordinate2D_Vec_dev_ptr);
        FracID = EleToFracID_ptr[EleID - 1];                              // from 0
        EdgeNO = make_uint3(EleID * 3 - 3, EleID * 3 - 2, EleID * 3 - 1); // from 0

        if (GGH == false)
        {
            printf("Warning: I cannot find which grid does the particle trajectory lies in! (I)\n");
            return;
        };

        goto StartOfCheckParticlePosition;
    };

    uint NumOfElesSharedEdge = EdgesSharedEle_DEV[CrossedGlobalEdgeNO].NumSharedEle;

    if (NumOfElesSharedEdge == 0)
    {
        printf("Warning: zero shared edge?\n");
        return;
    }

    if (NumOfElesSharedEdge == 1)
    {
        // we have to check if the edge is the target plane or not
        float RK[3][3];
        Frac_DEV[FracID].RoationMatrix(RK, 23);

        float3 IntersectionPnt_3D = make_float3(IntersectionPnt.x, IntersectionPnt.y, 0);
        IntersectionPnt_3D = cuDFNsys::ProductSquare3Float3(RK, IntersectionPnt_3D);
        IntersectionPnt_3D = make_float3(IntersectionPnt_3D.x + Frac_DEV[FracID].Center.x,
                                         IntersectionPnt_3D.y + Frac_DEV[FracID].Center.y,
                                         IntersectionPnt_3D.z + Frac_DEV[FracID].Center.z);

        float *tmp_coor_compo = &(IntersectionPnt_3D.x);

        if (abs(tmp_coor_compo[Dir_flow] - outletcoordinate) < 1e-3)
        {
            P_DEV[i].IfReachOutletPlane = true;
            P_DEV[i].Position2D = IntersectionPnt;

            return;
        };

        // if it is an inlet plane
        if (tmp_coor_compo[Dir_flow] > (-1.0f * outletcoordinate))
        {
            printf("Warning: a particle goes above the model!\n");
            return;
        }

        // if it is a non-flux plane
        if (tmp_coor_compo[Dir_flow] > outletcoordinate && tmp_coor_compo[Dir_flow] < (-1.0f * outletcoordinate))
        {
            float2 RefletedPosition = cuDFNsys::ParticleReflection(TargPos, Vertex_Triangle[(uint)ResultParticle.z], Vertex_Triangle[((uint)ResultParticle.z + 1) % 3]);
            TargPos = RefletedPosition;
            InitPos = IntersectionPnt;
            goto StartOfCheckParticlePosition;
        };
    }
    else if (NumOfElesSharedEdge == 2)
    {
        EleID = (EdgesSharedEle_DEV[CrossedGlobalEdgeNO].EleID[0] == EleID ? EdgesSharedEle_DEV[CrossedGlobalEdgeNO].EleID[1] : EdgesSharedEle_DEV[CrossedGlobalEdgeNO].EleID[0]);
        P_DEV[i].ElementID = EleID;
        InitPos = IntersectionPnt;

        FracID = EleToFracID_ptr[EleID - 1];                              // from 0
        EdgeNO = make_uint3(EleID * 3 - 3, EleID * 3 - 2, EleID * 3 - 1); // from 0

        // printf("StartOfCheckParticlePosition II, particle %d, eleID changed to %d\n\n", i + 1, EleID);

        if (EleID == Record_1stEleID && count_loop >= 5)
        {
            //printf("Warning: the detection of particle position may get stucked in am endless loop!\n");

            bool GGH = cuDFNsys::IdentifyParticleElementAllElements(InitPos, TargPos, EleID, numElements, EleToFracID_ptr,
                                                                    FracID, NeighborEleOfOneEle_dev_ptr[P_DEV[i].ElementID - 1], Coordinate2D_Vec_dev_ptr);
            FracID = EleToFracID_ptr[EleID - 1];                              // from 0
            EdgeNO = make_uint3(EleID * 3 - 3, EleID * 3 - 2, EleID * 3 - 1); // from 0

            // float2 Vertex_Triangle_PPP[3];
            // Vertex_Triangle_PPP[0] = make_float2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[0], Coordinate2D_Vec_dev_ptr[EleID - 1].y[0]);
            // Vertex_Triangle_PPP[1] = make_float2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[1], Coordinate2D_Vec_dev_ptr[EleID - 1].y[1]);
            // Vertex_Triangle_PPP[2] = make_float2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[2], Coordinate2D_Vec_dev_ptr[EleID - 1].y[2]);
            // printf("\nEleID == Record_1stEleID:\nParticleID %d, eleID: %d; Record_1stEleID: %d\nPntTrajectory:\n%f, %f\n%f, %f\nTri:\n%f, %f\n%f, %f\n%f, %f\nIntersection:\n%f, %f\nEdgeNO: %f\n\n",
            //        i + 1, EleID, Record_1stEleID,
            //        InitPos.x, InitPos.y, TargPos.x, TargPos.y,
            //        Vertex_Triangle_PPP[0].x, Vertex_Triangle_PPP[0].y,
            //        Vertex_Triangle_PPP[1].x, Vertex_Triangle_PPP[1].y,
            //        Vertex_Triangle_PPP[2].x, Vertex_Triangle_PPP[2].y,
            //        IntersectionPnt.x, IntersectionPnt.y, ResultParticle.z);

            if (GGH == false)
            {
                printf("Warning: I cannot find which grid does the particle trajectory lies in! (II)\n");
                return;
            }
        }

        goto StartOfCheckParticlePosition;
    }
    else
    {
        printf("Warning: more than two shared edge?\n");
        return;
        //
        uint PreviousFracID = FracID;

        // determine where to go
        cuDFNsys::WhichElementToGo(EleID,
                                   NumOfElesSharedEdge,
                                   EdgesSharedEle_DEV[CrossedGlobalEdgeNO].EleID,
                                   EleToFracID_ptr,
                                   Coordinate2D_Vec_dev_ptr,
                                   velocity_ptr,
                                   curand_uniform(&state),
                                   EleID,
                                   FracID);
        EdgeNO = make_uint3(EleID * 3 - 3, EleID * 3 - 2, EleID * 3 - 1); // from 0

        // change initial position
        float3 InitiPos3D = make_float3(IntersectionPnt.x, IntersectionPnt.y, 0.0f);
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

        // change target and intersection position
        // change target and intersection position
        float2 tmpDis = make_float2(TargPos.x - IntersectionPnt.x, TargPos.y - IntersectionPnt.y);
        float normDistance = sqrt(tmpDis.x * tmpDis.x + tmpDis.y + tmpDis.y);

        float3 Intersection3D = make_float3(IntersectionPnt.x, IntersectionPnt.y, 0),
               Target3D = make_float3(TargPos.x, TargPos.y, 0);

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
        float3 V_new = cuDFNsys::ProjectVToPlaneN(V_, Frac_DEV[FracID].NormalVec);
        V_new.x *= normDistance;
        V_new.y *= normDistance;
        V_new.z *= normDistance;

        ////
    };
};
}; // namespace cuDFNsys
