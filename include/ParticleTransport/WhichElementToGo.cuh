///////////////////////////////////////////////////////////////////
// NAME:              WhichElementToGo.cuh
//
// PURPOSE:           determine which element to go when the particle
//                    encounters an intersection
//
// FUNCTIONS/OBJECTS: WhichElementToGo
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../DataTypeSelector/DataTypeSelector.cuh"
#include "../GlobalDef/GlobalDef.cuh"
#include "../MHFEM/ReconstructVelocityGrid.cuh"
#include "../Mesh/EleCoor.cuh"

namespace cuDFNsys
{
template <typename T>
__host__ __device__ void WhichElementToGo(uint currentEleID,
                                          uint NumSharedEle,
                                          T Dispersion_local,
                                          uint *EleID_vec,
                                          uint *LocalEdgeNo_vec,
                                          uint *EleToFracID_ptr,
                                          cuDFNsys::EleCoor<T> *Coordinate2D_Vec_dev_ptr,
                                          T *velocity_ptr,
                                          T rand_0_1,
                                          int &NextElementID,
                                          int &NextFracID,
                                          int &IndexInLocal,
                                          bool &ifAllsharedEdgeVelocityPositive)
{
    T TotalVeloc = 0;
    T veloc_vec[_NumOfSharedEleAtMost];
    uint i_prime = 0;
    uint eleID_vec[_NumOfSharedEleAtMost];

    uint IndexTrans_vec[_NumOfSharedEleAtMost];

    for (uint i = 0; i < NumSharedEle; ++i)
    {
        if (EleID_vec[i] == currentEleID)
            continue;

        uint EleID = EleID_vec[i];               // from 1
        uint LocalEdgeNO__ = LocalEdgeNo_vec[i]; // 0, 1 or 2
        //uint FracID = EleToFracID_ptr[EleID - 1];                               // from 0
        uint3 EdgeNO = make_uint3(EleID * 3 - 3, EleID * 3 - 2, EleID * 3 - 1); // from 0

        cuDFNsys::Vector2<T> CenterThisTriangle = cuDFNsys::MakeVector2((T)1.0f / (T)3.0f * (Coordinate2D_Vec_dev_ptr[EleID - 1].x[0] + Coordinate2D_Vec_dev_ptr[EleID - 1].x[1] + Coordinate2D_Vec_dev_ptr[EleID - 1].x[2]),
                                                                        (T)1.0f / (T)3.0f * (Coordinate2D_Vec_dev_ptr[EleID - 1].y[0] + Coordinate2D_Vec_dev_ptr[EleID - 1].y[1] + Coordinate2D_Vec_dev_ptr[EleID - 1].y[2]));

        cuDFNsys::Vector3<T> Veloc_triangle = cuDFNsys::MakeVector3(velocity_ptr[EdgeNO.x],
                                                                    velocity_ptr[EdgeNO.y],
                                                                    velocity_ptr[EdgeNO.z]);
        T *tmpVelocity = &(Veloc_triangle.x);
        T VelocitySharedEdgeSep = tmpVelocity[LocalEdgeNO__];

        // printf("velocity: %f, %f, %f;\n", Veloc_triangle.x, Veloc_triangle.y, Veloc_triangle.z);
        // printf("normal Velocity of Shared Edge: %f\n", VelocitySharedEdgeSep);
        // printf("elementID: %d\n\n", EleID);

        if (VelocitySharedEdgeSep > 0)
            veloc_vec[i_prime] = 0;
        else
        {
            cuDFNsys::Vector2<T> Vertex_Triangle_ForVelocity[3];
            Vertex_Triangle_ForVelocity[0] = cuDFNsys::MakeVector2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[0], Coordinate2D_Vec_dev_ptr[EleID - 1].y[0]);
            Vertex_Triangle_ForVelocity[1] = cuDFNsys::MakeVector2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[1], Coordinate2D_Vec_dev_ptr[EleID - 1].y[1]);
            Vertex_Triangle_ForVelocity[2] = cuDFNsys::MakeVector2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[2], Coordinate2D_Vec_dev_ptr[EleID - 1].y[2]);

            cuDFNsys::Vector2<T> Veloc_p = cuDFNsys::ReconstructVelocityGrid<T>(CenterThisTriangle, Vertex_Triangle_ForVelocity, Veloc_triangle);
            // printf("velocity center: %f, %f;\n", Veloc_p.x, Veloc_p.y);

            T norm_veloc = sqrt(Veloc_p.x * Veloc_p.x + Veloc_p.y * Veloc_p.y);

            veloc_vec[i_prime] = norm_veloc;

            TotalVeloc += norm_veloc;
        };

        eleID_vec[i_prime] = EleID;
        IndexTrans_vec[i_prime] = LocalEdgeNO__;
        i_prime++;
    }

    ifAllsharedEdgeVelocityPositive = false;

    //printf("Dispersion_local: %.40f, TotalVeloc: %.40f\n", Dispersion_local,TotalVeloc);
    if (Dispersion_local == 0 && // particle tracking
        TotalVeloc == 0)         // all element normal velocities are positive
    {
        ifAllsharedEdgeVelocityPositive = true;
        return;
    }

    for (uint i = 0; i < NumSharedEle - 1; ++i)
    {
        veloc_vec[i] /= TotalVeloc;
        if (i > 0)
            veloc_vec[i] += veloc_vec[i - 1];
    }

    if (Dispersion_local > 0)
    {
        for (uint i = 0; i < NumSharedEle - 1; ++i)
            veloc_vec[i] = 1.0 / (NumSharedEle - 1) + (i > 0 ? veloc_vec[i - 1] : 0);
    }

    // printf("element velocity weight: ");
    // for (uint i = 0; i < NumSharedEle - 1; ++i)
    //     printf("%f, ", veloc_vec[i]);
    // printf("\n\n");

    for (uint i = 0; i < NumSharedEle - 1; ++i)
        if (rand_0_1 < veloc_vec[i])
        {
            NextElementID = eleID_vec[i];
            NextFracID = EleToFracID_ptr[NextElementID - 1];
            IndexInLocal = IndexTrans_vec[i];
            break;
        }
};
template __host__ __device__ void WhichElementToGo<double>(uint currentEleID,
                                                           uint NumSharedEle,
                                                           double Dispersion_local,
                                                           uint *EleID_vec,
                                                           uint *LocalEdgeNo_vec,
                                                           uint *EleToFracID_ptr,
                                                           cuDFNsys::EleCoor<double> *Coordinate2D_Vec_dev_ptr,
                                                           double *velocity_ptr,
                                                           double rand_0_1,
                                                           int &NextElementID,
                                                           int &NextFracID,
                                                           int &IndexInLocal,
                                                           bool &ifAllsharedEdgeVelocityPositive);
template __host__ __device__ void WhichElementToGo<float>(uint currentEleID,
                                                          uint NumSharedEle,
                                                          float Dispersion_local,
                                                          uint *EleID_vec,
                                                          uint *LocalEdgeNo_vec,
                                                          uint *EleToFracID_ptr,
                                                          cuDFNsys::EleCoor<float> *Coordinate2D_Vec_dev_ptr,
                                                          float *velocity_ptr,
                                                          float rand_0_1,
                                                          int &NextElementID,
                                                          int &NextFracID,
                                                          int &IndexInLocal,
                                                          bool &ifAllsharedEdgeVelocityPositive);
}; // namespace cuDFNsys