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

#include "ParticleTransport/Roate2DPositionTo3D.cuh"

// ====================================================
// NAME:        Roate2DPositionTo3D
// DESCRIPTION: Rotate a 2D particle position to 3D
// AUTHOR:      Tingchang YIN
// DATE:        15/11/2022
// ====================================================
template <typename T>
__host__ __device__ cuDFNsys::Vector3<T> cuDFNsys::Roate2DPositionTo3D(cuDFNsys::Vector2<T> PositionP,
                                                                       cuDFNsys::Fracture<T> OneFrac)
{
    cuDFNsys::Vector3<T> P_3D = cuDFNsys::MakeVector3(PositionP.x, PositionP.y, (T)0.0);

    T RK_2[3][3];
    OneFrac.RoationMatrix(RK_2, 23);

    P_3D = cuDFNsys::ProductSquare3Vector3<T>(RK_2, P_3D);
    P_3D.x += OneFrac.Center.x;
    P_3D.y += OneFrac.Center.y;
    P_3D.z += OneFrac.Center.z;

    return P_3D;
}; // Roate2DPositionTo3D
template __host__ __device__ cuDFNsys::Vector3<double> cuDFNsys::Roate2DPositionTo3D<double>(cuDFNsys::Vector2<double> PositionP,
                                                                                             cuDFNsys::Fracture<double> OneFrac);
template __host__ __device__ cuDFNsys::Vector3<float> cuDFNsys::Roate2DPositionTo3D<float>(cuDFNsys::Vector2<float> PositionP,
                                                                                           cuDFNsys::Fracture<float> OneFrac);