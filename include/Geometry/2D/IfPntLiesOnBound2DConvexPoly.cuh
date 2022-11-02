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
// NAME:              IfPntLiesOnBound2DConvexPoly.cuh
//
// PURPOSE:           if a point lies on bound of 2D convex polygon 2D
//
// FUNCTIONS/OBJECTS: IfPntLiesOnBound2DConvexPoly
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../../DataTypeSelector/DataTypeSelector.cuh"
#include "../../GlobalDef/GlobalDef.cuh"
#include "../../MatrixManipulation/MatrixManipulation.cuh"
#include "DistancePnt2DSeg.cuh"

namespace cuDFNsys
{
//if a point lies on 2D bound
template <typename T>
__device__ __host__ bool IfPntLiesOnBound2DConvexPoly(cuDFNsys::Vector2<T> pnt,
                                                      cuDFNsys::Vector2<T> *verts,
                                                      int N,
                                                      T _tol_);

// if a point lies on 2D bound, also return the edge NO
template <typename T>
__device__ __host__ bool IfPntLiesOnBound2DConvexPolyReturnEdgeNO(cuDFNsys::Vector2<T> pnt,
                                                                  cuDFNsys::Vector2<T> *verts,
                                                                  int N,
                                                                  T _tol_,
                                                                  int *edgeNO);
}; // namespace cuDFNsys