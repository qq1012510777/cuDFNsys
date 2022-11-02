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
// NAME:              IdentifyEleFrac.cuh
//
// PURPOSE:           Identify fracture ID of elements
//
// FUNCTIONS/OBJECTS: IdentifyEleFrac
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../Fractures/Fracture.cuh"
#include "../Geometry/2D/IfPntInside2DConvexPoly.cuh"
#include "../Geometry/2D/IfPntLiesOnBound2DConvexPoly.cuh"
#include "../Geometry/3D/DistancePnt3DPlane.cuh"
#include "../Geometry/3D/Triangle3DArea.cuh"
#include "../Geometry/3D/Scale3DTriangle.cuh"
#include "../GlobalDef/GlobalDef.cuh"

namespace cuDFNsys
{
template <typename T>
__global__ void IdentifyEleFrac(uint3 *One_entity_one_ele_dev_ptr,
                                cuDFNsys::Vector3<T> *coordinate_3D_dev_ptr,
                                cuDFNsys::Fracture<T> *Frac_verts_device_ptr,
                                int *Elements_Frac_dev_ptr,
                                int entity_count,
                                int frac_count,
                                T _tol_);
}; // namespace cuDFNsys