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

#include "Geometry/2D/Triangle2DOrientation.cuh"

// ====================================================
// NAME:        Triangle2DOrientation
// DESCRIPTION: return orientation of a 2D triangle
//              true: clockwise
//              false: counterclockwise
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================
template <typename T>
__device__ __host__ bool cuDFNsys::Triangle2DOrientation(cuDFNsys::Vector2<T> a,
                                                         cuDFNsys::Vector2<T> b,
                                                         cuDFNsys::Vector2<T> c)
{
    T d = (b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y);

    if (d < 0)
        return true;
    else
        return false;
}; // Triangle2DOrientation
template __device__ __host__ bool cuDFNsys::Triangle2DOrientation<double>(cuDFNsys::Vector2<double> a,
                                                                          cuDFNsys::Vector2<double> b,
                                                                          cuDFNsys::Vector2<double> c);
template __device__ __host__ bool cuDFNsys::Triangle2DOrientation<float>(cuDFNsys::Vector2<float> a,
                                                                         cuDFNsys::Vector2<float> b,
                                                                         cuDFNsys::Vector2<float> c);