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

#include "Geometry/3D/DistancePnt3DPlane.cuh"

// ====================================================
// NAME:        DistancePnt3DPlane
// DESCRIPTION: Distance between a 3D point and a 3D plane
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================
template <typename T>
__device__ __host__ T cuDFNsys::DistancePnt3DPlane(cuDFNsys::Vector3<T> Plane[3], cuDFNsys::Vector3<T> pnt)
{
    T a, b, c, d;

    T x1 = Plane[0].x;
    T x2 = Plane[1].x;
    T x3 = Plane[2].x;

    T y1 = Plane[0].y;
    T y2 = Plane[1].y;
    T y3 = Plane[2].y;

    T z1 = Plane[0].z;
    T z2 = Plane[1].z;
    T z3 = Plane[2].z;

    T a1 = x2 - x1;
    T b1 = y2 - y1;
    T c1 = z2 - z1;
    T a2 = x3 - x1;
    T b2 = y3 - y1;
    T c2 = z3 - z1;

    a = b1 * c2 - b2 * c1;
    b = a2 * c1 - a1 * c2;
    c = a1 * b2 - b1 * a2;
    d = (-a * x1 - b * y1 - c * z1);

    d = fabs((a * pnt.x + b * pnt.y +
              c * pnt.z + d));

    T e = sqrt(a * a + b * b + c * c);

    return (d / e);
}; // DistancePnt3DPlane
template __device__ __host__ double cuDFNsys::DistancePnt3DPlane<double>(cuDFNsys::Vector3<double> Plane[3], cuDFNsys::Vector3<double> pnt);
template __device__ __host__ float cuDFNsys::DistancePnt3DPlane<float>(cuDFNsys::Vector3<float> Plane[3], cuDFNsys::Vector3<float> pnt);