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
// NAME:              IdentifyPercolationCluster.cuh
//
// PURPOSE:           Identify Percolation Cluster,
//
// FUNCTIONS/OBJECTS: IdentifyPercolationCluster
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////

#pragma once
#include "../DataTypeSelector/DataTypeSelector.cuh"
#include "../GlobalDef/GlobalDef.cuh"
#include "Fracture.cuh"

namespace cuDFNsys
{
template <typename T>
class IdentifyPercolationCluster
{
public:
    // constructor
    IdentifyPercolationCluster(const std::vector<std::vector<size_t>> &ListClusters,
                               const thrust::host_vector<cuDFNsys::Fracture<T>> &Frac_verts_host,
                               const int &dir,
                               std::vector<size_t> &Percolation_cluster);
};
}; // namespace cuDFNsys