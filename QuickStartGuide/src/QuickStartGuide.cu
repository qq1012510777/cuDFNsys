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
// NAME:        A quickstart example
// DESCRIPTION: Call cuDFNsys functions to do simulation.
//              read Manual to understand this code
// AUTHOR:      Tingchang YIN
// DATE:        19/10/2023
// ====================================================

#include "cuDFNsys.cuh"
#include <fstream>
#include <iostream>
#include <limits.h>
#include <unistd.h>

int main(int argc, char *argv[])
{
    int NumFractures = 500;
    double kappa = 0;
    double beta = 0.25, gamma = 5.5e-4;
    int ModeSizeDistri = 0;
    cuDFNsys::Vector4<double> ParaSizeDistri =
        cuDFNsys::MakeVector4(1.5,
                              1.,
                              15.,
                              0.);

    double DomainSize_X = 30;
    double3 DomainDimensionRatio = make_double3(1, 1, 2);
    int perco_dir = 2;
    thrust::host_vector<cuDFNsys::Fracture<double>> Frac_verts_host(NumFractures);
    thrust::device_vector<cuDFNsys::Fracture<double>> Frac_verts_device(NumFractures);
    cuDFNsys::Fracture<double> *Frac_verts_device_ptr;
    Frac_verts_device_ptr = thrust::raw_pointer_cast(Frac_verts_device.data());
    time_t t;
    time(&t);
    cuDFNsys::Fractures<double><<<NumFractures / 256 + 1, 256>>>(Frac_verts_device_ptr,
                                                                 (unsigned long)t,
                                                                 NumFractures,
                                                                 DomainSize_X,
                                                                 ModeSizeDistri,
                                                                 ParaSizeDistri,
                                                                 kappa,
                                                                 beta,
                                                                 gamma,
                                                                 DomainDimensionRatio);

    cudaDeviceSynchronize();
    Frac_verts_host = Frac_verts_device;

    std::map<pair<size_t, size_t>, pair<cuDFNsys::Vector3<double>, cuDFNsys::Vector3<double>>> Intersection_map;
    cuDFNsys::IdentifyIntersection<double> identifyInters{Frac_verts_host.size(),
                                                          Frac_verts_device_ptr,
                                                          true,
                                                          Intersection_map};
    std::vector<std::vector<size_t>> ListClusters;
    std::vector<size_t> Percolation_cluster;
    cuDFNsys::Graph<double> G{(size_t)NumFractures, Intersection_map};
    G.UseDFS(ListClusters);
    cuDFNsys::IdentifyPercolationCluster<double> IdentiClu{ListClusters,
                                                           Frac_verts_host,
                                                           perco_dir,
                                                           Percolation_cluster};

    if (Percolation_cluster.size() > 0)
    {
        std::vector<size_t> Fracs_percol;
        cuDFNsys::GetAllPercolatingFractures GetPer{Percolation_cluster,
                                                    ListClusters,
                                                    Fracs_percol};

        std::vector<pair<int, int>> IntersectionPair_percol;
        cuDFNsys::RemoveDeadEndFrac<double> RDEF{Fracs_percol,
                                                 IntersectionPair_percol,
                                                 (size_t)perco_dir,
                                                 Frac_verts_host,
                                                 Intersection_map,
                                                 false};
        cuDFNsys::Mesh<double> mesh{Frac_verts_host,
                                    IntersectionPair_percol,
                                    &Fracs_percol,
                                    1,
                                    3,
                                    perco_dir,
                                    DomainSize_X,
                                    DomainDimensionRatio};
        double mean_grid_area = mesh.MatlabPlot("DFN_mesh.h5",
                                                "DFN_mesh.m",
                                                Frac_verts_host,
                                                DomainSize_X,
                                                true,
                                                true,
                                                true,
                                                "DFN_mesh",
                                                DomainDimensionRatio);
        cuDFNsys::MHFEM<double> fem{mesh,
                                    Frac_verts_host,
                                    100,
                                    20,
                                    perco_dir,
                                    DomainSize_X,
                                    DomainDimensionRatio};
        double2 velocities = fem.MatlabPlot("MHFEM.h5",
                                            "MHFEM.m",
                                            Frac_verts_host,
                                            mesh,
                                            DomainSize_X,
                                            true,
                                            "MHFEM",
                                            DomainDimensionRatio);

        double meanV = velocities.x;
        double maxV = velocities.y;
        double meanTime = pow(mean_grid_area, 0.5) / maxV;

        int NumParticles = 10000;
        int NumTimeStep = 200;
        double Pe_number = 200;
        double LengthScale = 5.4772;
        double DiffusionMole = LengthScale / Pe_number * meanV;
        double Factor_mean_time_in_grid = 2;
        double DeltaT = meanTime / Factor_mean_time_in_grid; // delta T

        cuDFNsys::OutputObjectData<double> lk;
        lk.OutputFractures("FracturesForParticle.h5", Frac_verts_host, DomainSize_X);

        cuDFNsys::ParticleTransport<double> p{NumTimeStep,
                                              Frac_verts_host,
                                              mesh,
                                              fem,
                                              (uint)perco_dir,
                                              -0.5 * (&DomainDimensionRatio.x)[perco_dir] * DomainSize_X,
                                              NumParticles,
                                              DeltaT,
                                              DiffusionMole,
                                              "Particle_tracking",
                                              "Flux-weighted",
                                              "OutputAll",
                                              false,
                                              0,
                                              false,
                                              50,
                                              true,
                                              true,
                                              20,
                                              true};

        p.MatlabPlot("MHFEM.h5",
                     "ParticlesMovement.m",
                     mesh,
                     fem,
                     DomainSize_X,
                     DomainDimensionRatio,
                     true,
                     "ParticlesMovement");
    }
    return 0;
};