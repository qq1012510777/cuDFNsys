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

#include "Geometry/2D/OrientationThree2DPnts.cuh"
// ====================================================
// NAME:        OrientationThree2DPnts
// DESCRIPTION: return orientation of three 2D points
//                    0: collinear
//                    1: clockwise
//                    2: counterclockwise
// AUTHOR:      Tingchang YIN
// DATE:        07/05/2022
// ====================================================

template <typename T>
__device__ __host__ uint cuDFNsys::OrientationThree2DPnts(cuDFNsys::Vector2<T> p, cuDFNsys::Vector2<T> q, cuDFNsys::Vector2<T> r, T _tol_)
{
    T val = (q.y - p.y) * (r.x - q.x) -
            (q.x - p.x) * (r.y - q.y);

    if (abs(val) < _tol_)
        return 0; // collinear

    // clock or counterclock wise
    return (val > 0) ? 1 : 2;
}; // OrientationThree2DPnts
template __device__ __host__ uint cuDFNsys::OrientationThree2DPnts<double>(cuDFNsys::Vector2<double> p, cuDFNsys::Vector2<double> q, cuDFNsys::Vector2<double> r, double _tol_);
template __device__ __host__ uint cuDFNsys::OrientationThree2DPnts<float>(cuDFNsys::Vector2<float> p, cuDFNsys::Vector2<float> q, cuDFNsys::Vector2<float> r, float _tol_);