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
// NAME:              IfParticlePositionInNeighboringElement.cuh
//
// PURPOSE:           If the particle position is in a neighboring element
//
// FUNCTIONS/OBJECTS: IfParticlePositionInNeighboringElement
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../Geometry/Geometry.cuh"
#include "../Mesh/EleCoor.cuh"
#include "EdgeToEle.cuh"
#include "IfParticleOnBoundOfElement.cuh"
#include "IfParticlePositionNearOneVertexOfElement.cuh"
#include "NeighborEle.cuh"

namespace cuDFNsys
{
template <typename T>
__device__ __host__ bool IfParticlePositionInNeighboringElement(cuDFNsys::Vector2<T> Position_p,
                                                                uint &EleID,
                                                                uint *EleToFracID_ptr,
                                                                uint FracID,
                                                                cuDFNsys::NeighborEle NE,
                                                                cuDFNsys::EleCoor<T> *Coordinate2D_Vec_dev_ptr,
                                                                T _TOL_p);
}; // namespace cuDFNsys