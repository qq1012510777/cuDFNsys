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
// NAME:              AssembleOnGPUKernel.cuh
//
// PURPOSE:           kernel function to assemble matrix
//
// FUNCTIONS/OBJECTS: AssembleOnGPUKernel
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../DataTypeSelector/DataTypeSelector.cuh"
#include "../GlobalDef/GlobalDef.cuh"
#include "../Mesh/Mesh.cuh"
#include "StimaA.cuh"
#include "Triplet.cuh"

namespace cuDFNsys
{
    template <typename T>
    __global__ void AssembleOnGPUKernel(
        cuDFNsys::Triplet<T> *tri_dev, cuDFNsys::EleCoor<T> *coord_2D_dev,
        cuDFNsys::EleEdgeAttri *Edge_attri, T *Conduc_Frac_dev,
        //int *neuman_sep_dev,
        int NUM_sep_edges, int NUM_eles, int Dim, T P_in, T P_out,
        bool IfPeriodic, T mu_over_RhoGravity = (T)1., T const_q = 1e-5);
}; // namespace cuDFNsys