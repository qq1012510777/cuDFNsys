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

#include "Geometry/3D/Scale3DTriangle.cuh"

// ====================================================
// NAME:        Scale3DTriangle
// DESCRIPTION: Scale a 3D triangle
// AUTHOR:      Tingchang YIN
// DATE:        22/08/2022
// ====================================================
template <typename T>
__device__ __host__ void cuDFNsys::Scale3DTriangle(cuDFNsys::Vector3<T> Triangle[3],
                                                   T factor)
{
    cuDFNsys::Vector3<T> Center;
    Center.x = (Triangle[0].x + Triangle[1].x + Triangle[2].x) / 3.0;
    Center.y = (Triangle[0].y + Triangle[1].y + Triangle[2].y) / 3.0;
    Center.z = (Triangle[0].z + Triangle[1].z + Triangle[2].z) / 3.0;

    for (uint i = 0; i < 3; ++i)
    {
        cuDFNsys::Vector3<T> dir = cuDFNsys::MakeVector3(Triangle[i].x - Center.x,
                                                         Triangle[i].y - Center.y,
                                                         Triangle[i].z - Center.z);
        //T norm = sqrt(dir.x * dir.x + dir.y * dir.y + dir.z * dir.z);

        Triangle[i].x = Center.x + factor * dir.x;
        Triangle[i].y = Center.y + factor * dir.y;
        Triangle[i].z = Center.z + factor * dir.z;
    };
};
template __device__ __host__ void cuDFNsys::Scale3DTriangle<double>(cuDFNsys::Vector3<double> Triangle[3],
                                                                    double factor);
template __device__ __host__ void cuDFNsys::Scale3DTriangle<float>(cuDFNsys::Vector3<float> Triangle[3],
                                                                   float factor);