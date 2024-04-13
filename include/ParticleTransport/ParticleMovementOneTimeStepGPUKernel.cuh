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
// NAME:              ParticleMovementOneTimeStepGPUKernel.cuh
//
// PURPOSE:           GPU kernel function to move particles
//                    in one time step
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
#include "AngleBetweenTwoNeighboringTriangles.cuh"
#include "EdgeToEle.cuh"
#include "IdentifyParticleCrossesWhichEdge.cuh"
#include "IfParticleOnBoundOfElement.cuh"
#include "IfParticlePositionInNeighboringElement.cuh"
#include "IfParticlePositionNearOneVertexOfElement.cuh"
#include "Particle.cuh"
#include "ParticleReflection.cuh"
#include "RandomWalkerGeometryFunctions.cuh"
#include "Roate2DPositionTo3D.cuh"
#include "RotationRemainningTrajectoryBetweenTwo3DTriangles.cuh"
#include "WhichElementToGo.cuh"

namespace cuDFNsys
{
    template <typename T>
    __global__ void ParticleMovementOneTimeStepGPUKernel(
        unsigned long seed, T delta_T_, T Dispersion_local,
        cuDFNsys::Particle<T> *P_DEV, cuDFNsys::Fracture<T> *Frac_DEV,
        cuDFNsys::EdgeToEle *EdgesSharedEle_DEV,
        cuDFNsys::EleCoor<T> *Coordinate2D_Vec_dev_ptr,
        //cuDFNsys::NeighborEle *NeighborEleOfOneEle_dev_ptr,
        uint *EleToFracID_ptr, T *velocity_ptr, uint Dir_flow,
        T outletcoordinate, int count, int numElements, uint stepNO,
        uint *Particle_runtime_error_dev_pnt, uint NUMParticlesInTotal,
        bool If_completeMixing_fluxWeighted, bool If_periodic,
        int2 *CorrespondingEleLocalEdge_device_ptr);
}; // namespace cuDFNsys
