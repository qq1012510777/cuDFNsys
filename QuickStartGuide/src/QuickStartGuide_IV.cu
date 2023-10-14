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
// NAME:        A quickstart example to do particle tracking
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
    cuDFNsys::MHFEM<double> fem;
    lk.InputMHFEM("mhfem.h5", fem);

    std::vector<double> tempVariable = hdf5Class.ReadDataset<double>("DFN_mesh.h5", "N", "mean_grid_area");
    double mean_grid_area = tempVariable[0];
    tempVariable = hdf5Class.ReadDataset<double>("DFN_MHFEM.h5", "N", "maxV");
    double maxV = tempVariable[0];
    tempVariable = hdf5Class.ReadDataset<double>("DFN_MHFEM.h5", "N", "meanV");
    double meanV = tempVariable[0];

    double meanTime_crossGrid = pow(mean_grid_area, 0.5) / maxV;

    double LengthScale = DomainSize_X;

    int NumParticles = 10000;
    double Pe_number = 200;
    int NumTimeStep = 200; // number of time steps
    double DiffusionMole = LengthScale / Pe_number * meanV;

    double DeltaT = meanTime_crossGrid / 3;

    cout << "- NumParticles: " << NumParticles << endl;
    cout << "- NumTimeStep: " << NumTimeStep << endl;
    cout << "- Pe_number: " << Pe_number << endl;
    cout << "- LengthScale: " << LengthScale << endl;
    cout << "- DeltaT: " << DeltaT << endl;
    cout << "- DiffusionMole: " << DiffusionMole << endl;
    cout << "\n";

    lk_out.OutputFractures("FracturesForParticle.h5", Frac_host, DomainSize_X);
    // the above two command just in order to output fractures information to transform 2D particle data to 3D

    double istart = cuDFNsys::CPUSecond();
    cuDFNsys::ParticleTransport<double> p{NumTimeStep,                                                // number of time steps
                                          Frac_host,                                                  // fractures
                                          mesh,                                                       // mesh
                                          fem,                                                        // fem
                                          (uint)Perco_dir,                                            // percolation direction
                                          -0.5 * (&DomainDimensionRatio.x)[Perco_dir] * DomainSize_X, // the target plane, z = -0.5 * L; it could be changed if the percolation direction is not along the z axis
                                          NumParticles,                                               // number of particles
                                          DeltaT,                                                     // delta t
                                          DiffusionMole,                                              // molecular diffusion
                                          "Particle_tracking",                                        // use particle tracking algorithm
                                          "Flux-weighted",                                            // injection mode
                                          "OutputAll",                                                // output all information
                                          false,                                                      // if use CPU?
                                          0,                                                          // number of CPU processors
                                          false,                                                      // if record run time
                                          20000,
                                          true,
                                          true,
                                          0.4 * (&DomainDimensionRatio.x)[Perco_dir] * DomainSize_X,
                                          false};

    p.MatlabPlot("DFN_MHFEM.h5",        // h5 file of mhfem
                 "ParticlesMovement.m", // matlab script
                 mesh,                  // mesh result
                 fem,                   // mhfem object
                 DomainSize_X,          // domain size
                 DomainDimensionRatio,  // ratio of dimensions of the domain
                 false,                 // Visualize it with python
                 "N");                  // the name of python script without suffix

    double ielaps = cuDFNsys::CPUSecond() - istart;
    cout << "Running time of particle tracking: " << ielaps << " sec\n";

    return 0;
};