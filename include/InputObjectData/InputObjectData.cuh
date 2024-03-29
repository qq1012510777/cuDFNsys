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
// NAME:              InputObjectData.cuh
//
// PURPOSE:           Input data of fractures, mesh, and so on
//
// FUNCTIONS/OBJECTS: InputObjectData
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../DataTypeSelector/DataTypeSelector.cuh"
#include "../Fractures/Fracture.cuh"
#include "../GlobalDef/GlobalDef.cuh"
#include "../HDF5API/HDF5API.cuh"
#include "../Mesh/Mesh.cuh"
#include "../MHFEM/MHFEM.cuh"

namespace cuDFNsys
{
template <typename T>
class InputObjectData
{
public:
    // constructor
    InputObjectData(){};

public:
    // Input fractures' data
    void InputFractures(const string &filename,
                        thrust::host_vector<cuDFNsys::Fracture<T>> &Frac_verts_host,
                        T &L,
                        double3 &DomainDimensionRatio);

    // input mesh data
    void InputMesh(const string &filename,
                   cuDFNsys::Mesh<T> &mesh,
                   std::vector<size_t> *Fracs_percol);

    // input mhfem data
    void InputMHFEM(const string &filename,
                    cuDFNsys::MHFEM<T> &MHFEM);
};
}; // namespace cuDFNsys