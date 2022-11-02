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
// NAME:              Intersection3DPolyXYPlane.cuh
//
// PURPOSE:           Identify intersection between
//                    a 3D Polygon and the XY plane
//
// FUNCTIONS/OBJECTS: Intersection3DPolyXYPlane
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../../DataTypeSelector/DataTypeSelector.cuh"
#include "../../GlobalDef/GlobalDef.cuh"
#include "../../MatrixManipulation/MatrixManipulation.cuh"
#include "Intersection3DSegXYPlane.cuh"

namespace cuDFNsys
{
// Identify intersection between a 3D Polygon and the XY plane
template <typename T>
__device__ __host__ bool Intersection3DPolyXYPlane(cuDFNsys::Vector3<T> *Poly,
                                                   int NUM_vert,
                                                   cuDFNsys::Vector3<T> *Intersection,
                                                   T _TOL_);
}; // namespace cuDFNsys