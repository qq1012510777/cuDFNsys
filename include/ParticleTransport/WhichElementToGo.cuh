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
                                          cuDFNsys::Fracture<T> *Frac_DEV,
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

        uint EleID = EleID_vec[i]; // from 1

        /// if (currentEleID == 13510)
        ///     printf("%d, ", EleID);

        uint LocalEdgeNO__ = LocalEdgeNo_vec[i];  // 0, 1 or 2
        uint FracID = EleToFracID_ptr[EleID - 1]; // from 0

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
            T conductivity_k = Frac_DEV[FracID].Conductivity;
            T b_aperture = pow(conductivity_k * 12.0, 1.0 / 3.0);
            Veloc_p.x /= b_aperture, Veloc_p.y /= b_aperture;

            // printf("velocity (LT^{-1}) center: %.40f, %.40f;\n", Veloc_p.x, Veloc_p.y);

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
    //printf("element velocity weight 0: ");
    //for (uint i = 0; i < NumSharedEle - 1; ++i)
    //    printf("%.40f, ", veloc_vec[i] / TotalVeloc);
    //printf("\n");

    // 2022-10-27 added the condition
    T velocity_error = (T)0.05;
    T aty = 0;
    for (uint i = 0; i < NumSharedEle - 1; ++i)
    {
        T avs = veloc_vec[i] / TotalVeloc;

        if (avs < velocity_error)
        {          
            aty += veloc_vec[i];
            veloc_vec[i] = 0;
        }
    }
    TotalVeloc -= aty;

    for (uint i = 0; i < NumSharedEle - 1; ++i)
    {
        veloc_vec[i] /= TotalVeloc;
        if (i > 0)
            veloc_vec[i] += veloc_vec[i - 1];
    }

    //printf("element velocity weight 1: ");
    //for (uint i = 0; i < NumSharedEle - 1; ++i)
    //    printf("%.40f, ", veloc_vec[i]);
    //printf("\n");

    // 2022-10-27 commented
    /// if (Dispersion_local == 0)
    /// {
    ///     T velocity_error = (T)0.05;
    ///     for (uint i = 0; i < NumSharedEle - 2; ++i)
    ///     {
    ///         if (veloc_vec[i] - (i == 0 ? 0 : veloc_vec[i - 1]) < velocity_error)
    ///             veloc_vec[i] = (i == 0 ? 0 : veloc_vec[i - 1]);
    ///         if (veloc_vec[i] > 1.0 - velocity_error)
    ///             veloc_vec[i] = 1.0;
    ///     }
    /// };

    //if (currentEleID == 13510)
    //{
    //for (uint i = 0; i < NumSharedEle - 1; ++i)
    //    printf("%.40f, ", veloc_vec[i]);
    //printf("\n");
    //}

    if (Dispersion_local > 0)
    {
        for (uint i = 0; i < NumSharedEle - 1; ++i)
            veloc_vec[i] = 1.0 / (NumSharedEle - 1) + (i > 0 ? veloc_vec[i - 1] : 0);
    }

    //printf("element velocity weight 2: ");
    //for (uint i = 0; i < NumSharedEle - 1; ++i)
    //   printf("%.40f, ", veloc_vec[i]);
    //printf("\n\n");

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
                                                           cuDFNsys::Fracture<double> *Frac_DEV,
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
                                                          cuDFNsys::Fracture<float> *Frac_DEV,
                                                          cuDFNsys::EleCoor<float> *Coordinate2D_Vec_dev_ptr,
                                                          float *velocity_ptr,
                                                          float rand_0_1,
                                                          int &NextElementID,
                                                          int &NextFracID,
                                                          int &IndexInLocal,
                                                          bool &ifAllsharedEdgeVelocityPositive);
}; // namespace cuDFNsys