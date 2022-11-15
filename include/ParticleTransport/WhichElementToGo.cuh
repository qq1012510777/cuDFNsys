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
// NAME:              WhichElementToGo.cuh
//
// PURPOSE:           determine which element to go when the particle
//                    encounters an intersection
//
// FUNCTIONS/OBJECTS: WhichElementToGo
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../DataTypeSelector/DataTypeSelector.cuh"
#include "../GlobalDef/GlobalDef.cuh"
#include "../MHFEM/MHFEM.cuh"
#include "../MHFEM/ReconstructVelocityGrid.cuh"
#include "../Mesh/EleCoor.cuh"

namespace cuDFNsys
{
template <typename T>
__host__ __device__ void WhichElementToGo(uint currentEleID,
                                          uint NumSharedEle,
                                          T Dispersion_local,
                                          uint *EleID_vec,
                                          uint *LocalEdgeNo_vec,
                                          uint *EleToFracID_ptr,
                                          cuDFNsys::Fracture<T> *Frac_DEV,
                                          cuDFNsys::EleCoor<T> *Coordinate2D_Vec_dev_ptr,
                                          T *velocity_ptr,
                                          T rand_0_1,
                                          int &NextElementID,
                                          int &NextFracID,
                                          int &IndexInLocal,
                                          bool &ifAllsharedEdgeVelocityPositive);
}; // namespace cuDFNsys