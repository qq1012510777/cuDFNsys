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
// NAME:              Mesh.cuh
//
// PURPOSE:           mesh of DFN
//
// FUNCTIONS/OBJECTS: Mesh
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../CPUSecond/CPUSecond.cuh"
#include "../DataTypeSelector/DataTypeSelector.cuh"
#include "../Exceptions/Exceptions.cuh"
#include "../Fractures/Fracture.cuh"
#include "../Geometry/Geometry.cuh"
#include "../GlobalDef/GlobalDef.cuh"
#include "../HDF5API/HDF5API.cuh"
//#include "../MatlabAPI/MatlabAPI.cuh"
#include "EleCoor.cuh"
#include "EleEdgeAttri.cuh"
#include "GetLocalCoordiates.cuh"
#include "IdentifyEleFrac.cuh"
#include "UMapEdge.cuh"
#include "gmsh.h"
#include <cstring>
#include <fstream>
#include <sstream>
#include <unordered_map>

namespace cuDFNsys
{
template <typename T>
class Mesh
{
public:
    // if mesh generates successfully or not?
    bool MeshSuccess = true;
    // Fracture ID of each polygon
    std::vector<size_t> *FracID;

    // 3D coordinates
    thrust::host_vector<cuDFNsys::Vector3<T>> Coordinate3D;
    // 2D elements of each fracture
    thrust::host_vector<thrust::host_vector<uint3>> Element2D;
    // 3D elements of whole DFN
    thrust::host_vector<uint3> Element3D;
    // 2D coordinates of elements
    thrust::host_vector<cuDFNsys::EleCoor<T>> Coordinate2D;
    // Frac Tag of each element
    thrust::host_vector<uint> ElementFracTag; // from 0

    // edge attributes
    thrust::host_vector<cuDFNsys::EleEdgeAttri> EdgeAttri;
    // length values of inlet edges
    thrust::host_vector<cuDFNsys::Vector2<T>> InletEdgeNOLen;
    // length values of outlet edges
    thrust::host_vector<cuDFNsys::Vector2<T>> OutletEdgeNOLen;

    // number of interior edges
    uint NumInteriorEdges;
    // number of inlet edges
    uint NumInletEdges;
    // number of outlet edges
    uint NumOutletEdges;
    // number of neaumann edges
    uint NumNeumannEdges;

public:
    // checked percolation direction
    int Dir = 2;

public:
    // constructor
    Mesh(){};

    // constructor
    Mesh(const thrust::host_vector<cuDFNsys::Fracture<T>> &Fracs,
         const std::vector<pair<int, int>> &IntersectionPair_percol,
         std::vector<size_t> *Fracs_percol,
         const T &min_ele_edge,
         const T &max_ele_edge,
         const int &dir_,
         const T &L);

    // plot mesh
    void MatlabPlot(const string &mat_key,
                    const string &command_key,
                    thrust::host_vector<cuDFNsys::Fracture<T>> Fracs,
                    const T &L,
                    const bool &if_check_2D_coordinates,
                    const bool &if_check_edge_Numbering,
                    bool if_python_visualization = false,
                    string PythonName_Without_suffix = "DFN_mesh_py");

private:
    // get coordinates of mesh
    void GetCoordinates();
    // get elements of mesh
    void GetElements(const thrust::host_vector<cuDFNsys::Fracture<T>> &Fracs_s);
    // numbering edges of elements
    void NumberingEdges(const T L);

private:
    // get elements in each 3D surface entity
    void GetEntitiesElements(thrust::host_vector<thrust::host_vector<uint3>> &elementEntities_2D,
                             thrust::host_vector<uint> &Largest_ele);
    // generate sparse matrix to check edge attributes
    UMapEdge SparseMatEdgeAttri(uint i, bool if_change_ori);
    // get element ID
    size_t GetElementID(size_t i, size_t j);
    // check if a edge is dirchlet edge
    pair<bool, string> IfTwoEndsDirchlet(const size_t node1,
                                         const size_t node2,
                                         const T L);
};
}; // namespace cuDFNsys
