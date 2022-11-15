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

#include "ParticleTransport/IfParticlePositionNearOneVertexOfElement.cuh"

// ====================================================
// NAME:        IfParticlePositionNearOneVertexOfElement
// DESCRIPTION: If the particle is near a vertex of element
// AUTHOR:      Tingchang YIN
// DATE:        15/11/2022
// ====================================================
template <typename T>
__device__ __host__ bool cuDFNsys::IfParticlePositionNearOneVertexOfElement(cuDFNsys::Vector2<T> Position_p, cuDFNsys::Vector2<T> Tri[3], T _TOL_p)
{
    for (uint i = 0; i < 3; ++i)
    {
        cuDFNsys::Vector2<T> V;
        V.x = Position_p.x - Tri[i].x;
        V.y = Position_p.y - Tri[i].y;

        T norm_ = sqrt(V.x * V.x + V.y * V.y);

        if (norm_ < _TOL_p)
            return true;
    }
    return false;
};
template __device__ __host__ bool cuDFNsys::IfParticlePositionNearOneVertexOfElement<double>(cuDFNsys::Vector2<double> Position_p, cuDFNsys::Vector2<double> Tri[3], double _TOL_p);
template __device__ __host__ bool cuDFNsys::IfParticlePositionNearOneVertexOfElement<float>(cuDFNsys::Vector2<float> Position_p, cuDFNsys::Vector2<float> Tri[3], float _TOL_p);