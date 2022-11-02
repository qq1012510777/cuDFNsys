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

#include "Geometry/2D/Scale2DSegment.cuh"

// ====================================================
// NAME:        Scale2DSegment
// DESCRIPTION: scale a 2D line segment along the center
// AUTHOR:      Tingchang YIN
// DATE:        09/05/2022
// ====================================================
template <typename T>
__host__ __device__ void cuDFNsys::Scale2DSegment(cuDFNsys::Vector2<T> *Segment, T scaleF)
{
    cuDFNsys::Vector2<T> center_;
    center_.x = (T)0.5 * (Segment[0].x + Segment[1].x);
    center_.y = (T)0.5 * (Segment[0].y + Segment[1].y);

    Segment[0].x = center_.x + (Segment[0].x - center_.x) * scaleF;
    Segment[0].y = center_.y + (Segment[0].y - center_.y) * scaleF;

    Segment[1].x = center_.x + (Segment[1].x - center_.x) * scaleF;
    Segment[1].y = center_.y + (Segment[1].y - center_.y) * scaleF;
}; // Scale2DSegment
template __host__ __device__ void cuDFNsys::Scale2DSegment<double>(cuDFNsys::Vector2<double> *Segment, double scaleF);
template __host__ __device__ void cuDFNsys::Scale2DSegment<float>(cuDFNsys::Vector2<float> *Segment, float scaleF);