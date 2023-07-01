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
// NAME:              RandomWalkerGeometryFunctions.cuh
//
// PURPOSE:           some geometry-related functions for random walkers
//
// FUNCTIONS/OBJECTS: ParticleMovementOneTimeStepGPUKernel
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../DataTypeSelector/DataTypeSelector.cuh"
#include "../Fractures/Fracture.cuh"
#include "../Geometry/Geometry.cuh"
#include "../GlobalDef/GlobalDef.cuh"
#include "../MHFEM/ReconstructVelocityGrid.cuh"
#include "../MatrixManipulation/MatrixManipulation.cuh"
#include "../Mesh/EleCoor.cuh"
#include "../RandomFunction/RandomFunction.cuh"

// ramdom walker geometry functions
namespace cuDFNsys
{
template <typename T>
__host__ __device__ uint2 IfRandomWalkLieOnBorderTriangle(cuDFNsys::Vector2<T> p,
                                                          cuDFNsys::Vector2<T> Triangle[3],
                                                          T Tol);

template <typename T>
__host__ __device__ int Sgn(const T &x);

template <typename T>
__host__ __device__ bool WhichEdgeDoesTrajectoryIntersect(cuDFNsys::Vector2<T> pp[2],
                                                          cuDFNsys::Vector2<T> Triangle[3],
                                                          uint Result[4]);

template <typename T>
__host__ __device__ cuDFNsys::Vector3<T> IntersectionLineLine2D(cuDFNsys::Vector2<T> pp[2],
                                                                cuDFNsys::Vector2<T> cc[2]);

}; // namespace cuDFNsys
