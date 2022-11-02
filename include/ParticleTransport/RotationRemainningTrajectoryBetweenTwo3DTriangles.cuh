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
// NAME:              RotationRemainningTrajectoryBetweenTwo3DTriangles.cuh
//
// PURPOSE:           rotation remanning trajectory between two adjacent elements
//                    after the particle goes through a fracture trace
//
// FUNCTIONS/OBJECTS: RotationRemainningTrajectoryBetweenTwo3DTriangles
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../DataTypeSelector/DataTypeSelector.cuh"
#include "../GlobalDef/GlobalDef.cuh"
#include "../MatrixManipulation/MatrixManipulation.cuh"

namespace cuDFNsys
{
template <typename T>
__host__ __device__ cuDFNsys::Vector3<T> RotationRemainningTrajectoryBetweenTwo3DTriangles(const cuDFNsys::Vector3<T> preEle[3],
                                                                                           const T &theta_angle)
{
    cuDFNsys::Vector3<T> pntA = preEle[0],
                         pntB = preEle[1],
                         pntC = preEle[2];

    cuDFNsys::Vector3<T> BC = cuDFNsys::MakeVector3(pntC.x - pntB.x, pntC.y - pntB.y, pntC.z - pntB.z);
    T norm_V = sqrt(BC.x * BC.x + BC.y * BC.y + BC.z * BC.z);
    BC.x /= norm_V;
    BC.y /= norm_V;
    BC.z /= norm_V;

    cuDFNsys::Vector3<T> Pnt = cuDFNsys::MakeVector3(pntA.x - pntB.x, pntA.y - pntB.y, pntA.z - pntB.z);

    T sin_theta = sin(theta_angle);
    T cos_theta = cos(theta_angle);

    T RK[3][3] = {cos_theta + BC.x * BC.x * (1 - cos_theta), BC.x * BC.y * (1 - cos_theta) - BC.z * sin_theta, BC.x * BC.z * (1 - cos_theta) + BC.y * sin_theta,
                  BC.x * BC.y * (1 - cos_theta) + BC.z * sin_theta, cos_theta + BC.y * BC.y * (1 - cos_theta), BC.y * BC.z * (1 - cos_theta) - BC.x * sin_theta,
                  BC.x * BC.z * (1 - cos_theta) - BC.y * sin_theta, BC.y * BC.z * (1 - cos_theta) + BC.x * sin_theta, cos_theta + BC.z * BC.z * (1 - cos_theta)};

    cuDFNsys::Vector3<T> new_Pnt;
    new_Pnt = cuDFNsys::ProductSquare3Vector3<T>(RK, Pnt);
    new_Pnt.x += pntB.x;
    new_Pnt.y += pntB.y;
    new_Pnt.z += pntB.z;

    return new_Pnt;
};
template __host__ __device__ cuDFNsys::Vector3<double> RotationRemainningTrajectoryBetweenTwo3DTriangles<double>(const cuDFNsys::Vector3<double> preEle[3],
                                                                                                                 const double &theta_angle);
template __host__ __device__ cuDFNsys::Vector3<float> RotationRemainningTrajectoryBetweenTwo3DTriangles<float>(const cuDFNsys::Vector3<float> preEle[3],
                                                                                                               const float &theta_angle);
}; // namespace cuDFNsys