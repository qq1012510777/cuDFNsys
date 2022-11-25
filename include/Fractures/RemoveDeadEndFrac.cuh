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
// NAME:              RemoveDeadEndFrac.cuh
//
// PURPOSE:           Remove dead end fractures
//                    that do not conduct flow
//
// FUNCTIONS/OBJECTS: RemoveDeadEndFrac
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../DataTypeSelector/DataTypeSelector.cuh"
#include "../GlobalDef/GlobalDef.cuh"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Fracture.cuh"
#include <vector>

using namespace std;
using namespace Eigen;

namespace cuDFNsys
{
template <typename T>
class RemoveDeadEndFrac
{
public:
    RemoveDeadEndFrac(std::vector<size_t> &One_cluster,
                      std::vector<pair<int, int>> &Intersection_pair,
                      const size_t &dir,
                      thrust::host_vector<cuDFNsys::Fracture<T>> &Fracs,
                      std::map<pair<size_t, size_t>, pair<cuDFNsys::Vector3<T>, cuDFNsys::Vector3<T>>> Intersection_map,
                      bool IfRemoveDeadEnds = true);
};
}; // namespace cuDFNsys