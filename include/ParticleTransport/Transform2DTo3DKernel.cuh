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
// NAME:              Transform2DTo3DKernel.cuh
//
// PURPOSE:           Transform 2D particle positions to 3D
//
// FUNCTIONS/OBJECTS: Transform2DTo3DKernel
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "EdgeToEle.cuh"
#include "NeighborEle.cuh"
#include "OutputObjectData/OutputObjectData.cuh"
#include "Particle.cuh"
#include "ParticleMovementOneTimeStepCPU.cuh"
#include "ParticleMovementOneTimeStepGPUKernel.cuh"
#include "PredicateNumOfReachedOutletParticles.cuh"
#include "ToStringWithWidth/ToStringWithWidth.cuh"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

namespace cuDFNsys
{
template <typename T>
__global__ void Transform2DTo3DKernel(cuDFNsys::Fracture<T> *Frac_verts_device_ptr,
                                      T *Position3D_dev_ptr,
                                      cuDFNsys::Particle<T> *temp2Dpos_dev_ptr,
                                      uint *ElementFracTag_cuda_devptr,
                                      //uint *EleTag_device_ptr,
                                      uint count, T outletcoordinate, uint Dir);
}; // namespace cuDFNsys