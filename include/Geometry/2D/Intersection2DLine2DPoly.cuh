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
// NAME:              Intersection2DLine2DPoly.cuh
//
// PURPOSE:           Intersection between
//                    2D Line (infinite) and 2D Polygon
//
// FUNCTIONS/OBJECTS: Intersection2DLine2DPoly
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////

#pragma once
#include "../../DataTypeSelector/DataTypeSelector.cuh"
#include "../../GlobalDef/GlobalDef.cuh"
#include "../../MatrixManipulation/MatrixManipulation.cuh"
#include "Intersection2DLine2DSeg.cuh"

namespace cuDFNsys
{
// Intersection between 2D Line (infinite) and 2D Polygon
template <typename T>
__device__ __host__ bool Intersection2DLine2DPoly(cuDFNsys::Vector2<T> *Poly2D,
                                                  int NUM_verts,
                                                  cuDFNsys::Vector2<T> *Line,
                                                  cuDFNsys::Vector2<T> *intersection_k,
                                                  T _TOL_);
}; // namespace cuDFNsys