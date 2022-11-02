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

#include "Geometry/2D/Triangle2DArea.cuh"

// ====================================================
// NAME:        Triangle2DArea
// DESCRIPTION: get 2D triangle area
// AUTHOR:      Tingchang YIN
// DATE:        19/04/2022
// ====================================================
template <typename T>
__device__ __host__ T cuDFNsys::Triangle2DArea(cuDFNsys::Vector2<T> Pnt1,
                                               cuDFNsys::Vector2<T> Pnt2,
                                               cuDFNsys::Vector2<T> Pnt3)
{
    cuDFNsys::Vector2<T> Pntsets[3] = {Pnt1, Pnt2, Pnt3};

    T area = 0;

    T L[3] = {0};

    for (int i = 0; i < 3; ++i)
    {
        cuDFNsys::Vector2<T> KK = cuDFNsys::MakeVector2(Pntsets[i].x - Pntsets[(i + 1) % 3].x,
                                                        Pntsets[i].y - Pntsets[(i + 1) % 3].y);

        L[i] = pow(KK.x * KK.x + KK.y * KK.y, 0.5);
    }

    T P = (L[0] + L[1] + L[2]) * 0.5;
    area = pow(P * (P - L[0]) * (P - L[1]) * (P - L[2]), 0.5);

    return area;
}; // Triangle2DArea
template __device__ __host__ double cuDFNsys::Triangle2DArea(cuDFNsys::Vector2<double> Pnt1,
                                                             cuDFNsys::Vector2<double> Pnt2,
                                                             cuDFNsys::Vector2<double> Pnt3);
template __device__ __host__ float cuDFNsys::Triangle2DArea(cuDFNsys::Vector2<float> Pnt1,
                                                            cuDFNsys::Vector2<float> Pnt2,
                                                            cuDFNsys::Vector2<float> Pnt3);
