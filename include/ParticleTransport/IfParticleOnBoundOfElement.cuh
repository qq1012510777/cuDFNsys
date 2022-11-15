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
// NAME:              IfParticleOnBoundOfElement.cuh
//
// PURPOSE:           If particle lies on bound of element
//
// FUNCTIONS/OBJECTS: IfParticleOnBoundOfElement
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../DataTypeSelector/DataTypeSelector.cuh"
#include "../Geometry/Geometry.cuh"
#include "../GlobalDef/GlobalDef.cuh"

namespace cuDFNsys
{
template <typename T>
__host__ __device__ bool IfParticleOnBoundOfElement(cuDFNsys::Vector2<T> PositionP,
                                                    cuDFNsys::Vector2<T> GridVertex[3],
                                                    int &EdgeNOLocal,
                                                    T _TOL_);
}; // namespace cuDFNsys