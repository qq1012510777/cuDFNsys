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
// NAME:              IdentifyParticleCrossesWhichEdge.cuh
//
// PURPOSE:           Identify which edge does the particle cross
//                      -1: no
//                       0 1 2
//
// FUNCTIONS/OBJECTS: IdentifyParticleCrossesWhichEdge
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../DataTypeSelector/DataTypeSelector.cuh"
#include "../Geometry/Geometry.cuh"
#include "../GlobalDef/GlobalDef.cuh"
#include "../MHFEM/ReconstructVelocityGrid.cuh"
#include "EdgeToEle.cuh"
#include "Particle.cuh"

namespace cuDFNsys
{
template <typename T>
__host__ __device__ cuDFNsys::Vector3<T> IdentifyParticleCrossesWhichEdge(cuDFNsys::Vector2<T> *P_Trajectory,
                                                                          cuDFNsys::Vector2<T> *Vertex_Triangle,
                                                                          T _TOL_,
                                                                          cuDFNsys::Vector2<T> CrossedGlobalEdge[10][2],
                                                                          int CountCrossedGlobalEdge,
                                                                          uint stepNO,
                                                                          uint particleNO);
}; // namespace cuDFNsys