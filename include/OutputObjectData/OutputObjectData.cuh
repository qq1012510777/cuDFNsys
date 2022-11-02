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
// NAME:              OutputObjectData.cuh
//
// PURPOSE:           Output data of fractures, mesh, and so on
//
// FUNCTIONS/OBJECTS: OutputObjectData
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../DataTypeSelector/DataTypeSelector.cuh"
#include "../Fractures/Fracture.cuh"
#include "../GlobalDef/GlobalDef.cuh"
#include "../HDF5API/HDF5API.cuh"
#include "../Mesh/Mesh.cuh"

namespace cuDFNsys
{
template <typename T>
class OutputObjectData
{
public:
    // constructor
    OutputObjectData(){};

public:
    // output Fractures in H5
    void OutputFractures(const string &filename_,
                         const thrust::host_vector<cuDFNsys::Fracture<T>> &Frac_verts_host,
                         const T &L);

    // output Mesh in h5
    void OutputMesh(const string &filename_,
                    cuDFNsys::Mesh<T> mesh,
                    const std::vector<size_t> &Fracs_percol);
};
}; // namespace cuDFNsys