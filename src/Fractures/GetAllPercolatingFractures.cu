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

#include "Fractures/GetAllPercolatingFractures.cuh"

// ====================================================
// NAME:        GetAllPercolatingFractures
// DESCRIPTION: Get all fractures in percolating clusters
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================
cuDFNsys::GetAllPercolatingFractures::GetAllPercolatingFractures(const std::vector<size_t> &Percolation_cluster,
                                                                 const std::vector<std::vector<size_t>> &ListClusters,
                                                                 std::vector<size_t> &Fracs_percol)
{
    size_t fracNUM = 0;

    for (size_t i = 0; i < Percolation_cluster.size(); ++i)
        fracNUM += ListClusters[Percolation_cluster[i]].size();

    Fracs_percol.clear();
    Fracs_percol.reserve(fracNUM);

    for (size_t i = 0; i < Percolation_cluster.size(); ++i)
        Fracs_percol.insert(Fracs_percol.end(),
                            ListClusters[Percolation_cluster[i]].begin(),
                            ListClusters[Percolation_cluster[i]].end());
}; // GetAllPercolatingFractures
