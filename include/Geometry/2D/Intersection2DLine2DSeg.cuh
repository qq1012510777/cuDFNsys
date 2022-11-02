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
// NAME:              Intersection2DLine2DSeg.cuh
//
// PURPOSE:           Identify intersection between
//                    2D line (infinite) and 2D segment
//
// FUNCTIONS/OBJECTS: Intersection2DLine2DSeg
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////

#pragma once
#include "../../DataTypeSelector/DataTypeSelector.cuh"
#include "../../GlobalDef/GlobalDef.cuh"
#include "../../MatrixManipulation/MatrixManipulation.cuh"

namespace cuDFNsys
{
// Identify intersection between 2D line (infinite) and 2D segment
template <typename T>
__device__ __host__ bool Intersection2DLine2DSeg(cuDFNsys::Vector2<T> *Line,
                                                 cuDFNsys::Vector2<T> *Seg,
                                                 int *sign_, // 1, pnt; 2, seg; 3, none
                                                 cuDFNsys::Vector2<T> *intersection,
                                                 T _TOL_);
}; // namespace cuDFNsys
