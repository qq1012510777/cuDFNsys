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

// ====================================================
// NAME:        A quickstart example to generate mesh in
//              the DFN
// DESCRIPTION: Call cuDFNsys functions to do simulation.
// AUTHOR:      Tingchang YIN
// DATE:        13/10/2023
// ====================================================

#include "cuDFNsys.cuh"
#include <fstream>
#include <iostream>
#include <limits.h>
#include <unistd.h>

// a function to load variables (ListClusters, Percolation_cluster, Intersection_map) from .h5 files
// the three variables are created in QuickStartGuide_I
// now they are also used in this .cu file
void LoadPercolationInformation(std::vector<std::vector<size_t>> &ListClusters,
                                std::vector<size_t> &Percolation_cluster,
                                std::map<pair<size_t, size_t>, pair<cuDFNsys::Vector3<double>, cuDFNsys::Vector3<double>>> &Intersection_map);

int main(int argc, char *argv[])
{

    int dev = 0;
    // No. 0 GPU card
    GPUErrCheck(cudaSetDevice(dev));
    // try to use the first GPU card and check its availability

    // in `QuickStartGuide_I.cu`, we have generate a DFN which is percolative
    // let's generate mesh in this DFN
    int Perco_dir = 2;

    double DomainSize_X;
    double3 DomainDimensionRatio;
    thrust::host_vector<cuDFNsys::Fracture<double>> Frac_host;

    cuDFNsys::InputObjectData<double> lk;
    // a C++ class to load objects' information, e.g., fractures, mesh, fem
    cuDFNsys::OutputObjectData<double> lk_out;
    // a C++ class to output objects' information (e.g., mesh, fractures, fem)

    lk.InputFractures("Fractures.h5", Frac_host, DomainSize_X, DomainDimensionRatio);

    std::vector<std::vector<size_t>> ListClusters; // will store all fracture clusters
    std::vector<size_t> Percolation_cluster;       // will store the No. of percolating cluster
    std::map<pair<size_t, size_t>, pair<cuDFNsys::Vector3<double>, cuDFNsys::Vector3<double>>> Intersection_map;

    // load ListClusters, Percolation_cluster and Intersection_map from h5 files
    LoadPercolationInformation(ListClusters, Percolation_cluster, Intersection_map);

    // now let's mesh
    std::vector<size_t> Fracs_percol;
    // will store ID / No of fractures in percolating cluster

    double istart = cuDFNsys::CPUSecond();
    cuDFNsys::GetAllPercolatingFractures GetPer{Percolation_cluster,
                                                ListClusters,
                                                Fracs_percol};
    // this function is simple, just collecting fractures in the percolating clusters

    std::vector<pair<int, int>> IntersectionPair_percol;
    // will store the intersection pair of fractures in percolation clusters

    cuDFNsys::RemoveDeadEndFrac<double> RDEF{Fracs_percol,            // fractures' ID in percolating cluster
                                             IntersectionPair_percol, // intersection pair
                                             (size_t)Perco_dir,       // flow direction
                                             Frac_host,               // host vector of fractures
                                             Intersection_map,        // map of intersection
                                             false};                  // does not remove dead-end fractures; true = remove dead ends

    // the above function just sorts fractures and does not remove dead end fractures

    lk_out.OutputFractures("Fractures.h5", Frac_host, DomainSize_X, DomainDimensionRatio);

    double ielaps = cuDFNsys::CPUSecond() - istart;
    cout << "Running time of sorting fractures: " << ielaps << " sec\n";

    istart = cuDFNsys::CPUSecond();
    cuDFNsys::Mesh<double> mesh{Frac_host,               // host vector of fractures, after removing fractures of dead end
                                IntersectionPair_percol, // intersection pair
                                &Fracs_percol,           // fractures' ID in percolating cluster
                                1,                       // minimum grid size
                                3,                       // maximum grid size
                                Perco_dir,               // flow direction
                                DomainSize_X,            // domain size
                                DomainDimensionRatio};

    lk_out.OutputMesh("mesh.h5", mesh, Fracs_percol);

    // mesh finished, mean_grid_area is the mean area of finite elements
    double mean_grid_area = mesh.MatlabPlot("DFN_mesh.h5", // h5 file
                                            "DFN_mesh.m",  // name of matlab script
                                            Frac_host,     // fracture vector on the host side
                                            DomainSize_X,  // domain size
                                            true,          // if check 2D coordinates, because 3D fractures can be mapped to 2D plane
                                            true,          // if check the edge attributes? Neumann, Dirichlet?
                                            true,          // if I want to see mesh with Python?
                                            "DFN_mesh",    // name of python script without suffix
                                            DomainDimensionRatio);

    cuDFNsys::HDF5API hdf5Class;

    return 0;
};

//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//

// source of LoadPercolationInformation
void LoadPercolationInformation(std::vector<std::vector<size_t>> &ListClusters,
                                std::vector<size_t> &Percolation_cluster,
                                std::map<pair<size_t, size_t>, pair<cuDFNsys::Vector3<double>, cuDFNsys::Vector3<double>>> &Intersection_map)
{
    cuDFNsys::HDF5API hdf5Class;
    std::vector<int> NumClusters = hdf5Class.ReadDataset<int>("DFN_truncated.h5", "N", "NumClusters");
    ListClusters.resize(NumClusters[0]);
    std::vector<double> Temp_Variable = hdf5Class.ReadDataset<double>("DFN_truncated.h5", "N", "ListClusters");
    for (int i = 0; i < NumClusters[0]; ++i)
        ListClusters[i].reserve(Temp_Variable.size() / NumClusters[0]);
    for (int i = 0; i < Temp_Variable.size(); ++i)
        if (Temp_Variable[i] != -1)
            ListClusters[i % NumClusters[0]].push_back(Temp_Variable[i] - 1);
    for (int i = 0; i < NumClusters[0]; ++i)
        ListClusters[i].shrink_to_fit();
    Temp_Variable = hdf5Class.ReadDataset<double>("DFN_truncated.h5", "N", "PercolationClusters");
    Percolation_cluster.resize(Temp_Variable.size());
    for (int i = 0; i < Temp_Variable.size(); ++i)
        Percolation_cluster[i] = Temp_Variable[i] - 1;

    //------INTERSECTION PAIR---
    Temp_Variable = hdf5Class.ReadDataset<double>("DFN_truncated.h5", "N", "intersections");
    int NumIntersections = Temp_Variable.size() / 8;
    for (int i = 0; i < NumIntersections; ++i)
        Intersection_map.insert(std::make_pair(std::make_pair((uint)Temp_Variable[i + 6 * NumIntersections], (uint)Temp_Variable[i + 7 * NumIntersections]),
                                               std::make_pair(cuDFNsys::MakeVector3(Temp_Variable[i], Temp_Variable[i + NumIntersections], Temp_Variable[i + NumIntersections * 2]), cuDFNsys::MakeVector3(Temp_Variable[i + NumIntersections * 3], Temp_Variable[i + NumIntersections * 4], Temp_Variable[i + NumIntersections * 5]))));
};