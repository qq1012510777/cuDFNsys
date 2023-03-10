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
// NAME:              MatlabPlotDFN.cuh
//
// PURPOSE:           Plot DFN with m and mat file
//                    m is command file; mat is data file
//                    updated: this function can also output py and h5 file
//
// FUNCTIONS/OBJECTS: MatlabPlotDFN
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../DataTypeSelector/DataTypeSelector.cuh"
#include "../GlobalDef/GlobalDef.cuh"
#include "../HDF5API/HDF5API.cuh"
//#include "../MatlabAPI/MatlabAPI.cuh"
#include "Fracture.cuh"
#include <algorithm>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace cuDFNsys
{
template <typename T>
class MatlabPlotDFN
{
public:
    // constructor
    MatlabPlotDFN(string mat_key,                                                                                     // mat file name
                  string command_key,                                                                                 // m file name
                  thrust::host_vector<cuDFNsys::Fracture<T>> Frac_verts_host,                                         // Vector of Fracture
                  std::map<pair<size_t, size_t>, pair<cuDFNsys::Vector3<T>, cuDFNsys::Vector3<T>>> Intersection_host, // vector of Intersection
                  std::vector<std::vector<size_t>> ListClusters,                                                      // List of clusters: a cluster is a vector
                  std::vector<size_t> Percolation_cluster,                                                            // Percolation cluster contains the cluster vector ID
                  bool If_show_truncated_frac,                                                                        // if show truncated fractures
                  bool If_show_intersection,
                  bool If_show_cluster,
                  bool If_show_orientation,
                  T L,
                  int dir,
                  bool if_python_visualization = false,
                  string PythonName_Without_suffix = "DFN_py",
                  double3 DomainDimensionRatio_d = make_double3(1, 1, 1));
};
}; // namespace cuDFNsys