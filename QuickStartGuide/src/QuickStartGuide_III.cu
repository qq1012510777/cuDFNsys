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
// NAME:        A quickstart example to solve flow in DFNs
// DESCRIPTION: Call cuDFNsys functions to do simulation.
// AUTHOR:      Tingchang YIN
// DATE:        13/10/2023
// ====================================================

#include "cuDFNsys.cuh"
#include <fstream>
#include <iostream>
#include <limits.h>
#include <unistd.h>

int main(int argc, char *argv[])
{

    int dev = 0;
    // No. 0 GPU card
    GPUErrCheck(cudaSetDevice(dev));
    // try to use the first GPU card and check its availability

    int Perco_dir = 2;
    double DomainSize_X;
    double3 DomainDimensionRatio;
    thrust::host_vector<cuDFNsys::Fracture<double>> Frac_host;

    cuDFNsys::HDF5API hdf5Class;

    cuDFNsys::InputObjectData<double> lk;
    // a C++ class to load objects' information, e.g., fractures, mesh, fem
    cuDFNsys::OutputObjectData<double> lk_out;
    // a C++ class to output objects' information (e.g., mesh, fractures, fem)

    lk.InputFractures("Fractures.h5", Frac_host, DomainSize_X, DomainDimensionRatio);

    std::vector<uint> Fracs_percol_II = hdf5Class.ReadDataset<uint>("mesh.h5", "N", "Fracs_percol");
    std::vector<size_t> Fracs_percol(Fracs_percol_II.size());
    std::copy(Fracs_percol_II.begin(), Fracs_percol_II.end(), Fracs_percol.data());
    cuDFNsys::Mesh<double> mesh;
    lk.InputMesh("mesh.h5", mesh, &Fracs_percol);

    cuDFNsys::MHFEM<double> fem{mesh,         // mesh object
                                Frac_host,    // fractures
                                DomainSize_X, // hydraulic head of inlet = 100 m
                                0,            // hydraulic head of outlet = 20 m
                                Perco_dir,    // flow direction
                                DomainSize_X, // domain size
                                DomainDimensionRatio};

    cout << "Fluxes, inlet: " << fem.QIn << ", outlet: ";
    cout << fem.QOut << ", Permeability: ";
    cout << fem.Permeability << endl;
    if (fem.QError > 1 || isnan(fem.Permeability) == 1)
        throw cuDFNsys::ExceptionsIgnore("Found large error or isnan, the error: " + std::to_string(fem.QError) + ", the permeability: " + std::to_string(fem.Permeability) + "\n");

    lk_out.OutputMHFEM("mhfem.h5", fem);

    double2 TGH = fem.MatlabPlot("DFN_MHFEM.h5", // h5 file
                                 "DFN_MHFEM.m",  // matlab script to see mhfem result
                                 Frac_host,      // fractures
                                 mesh,           // mesh object
                                 DomainSize_X,   // domain size
                                 true,           // if use python to do visualization
                                 "DFN_MHFEM",    // name of python script, without suffix
                                 DomainDimensionRatio);

    return 0;
};