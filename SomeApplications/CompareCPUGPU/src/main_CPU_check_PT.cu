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
// NAME:        A test case
// DESCRIPTION: Call cuDFNsys functions to do simulation and test.
// AUTHOR:      Tingchang YIN
// DATE:        30/06/2022
// ====================================================

#include "cuDFNsys.cuh"
#include <bits/stdc++.h>
#include <iostream>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#ifdef USE_DOUBLES
typedef double _DataType_;
#else
typedef float _DataType_;
#endif

const int LoopTime = 5;

int main(int argc, char *argv[])
{

    try
    {

        double istart = cuDFNsys::CPUSecond();

        int dev = 0;
        GPUErrCheck(cudaSetDevice(dev));

        int DSIZE = 4;
        _DataType_ L = 30;
        _DataType_ minGrid = 2.0;
        _DataType_ maxGrid = 5.0;

        int perco_dir = 2;

        cout << "preparing" << endl;

        thrust::host_vector<cuDFNsys::Fracture<_DataType_>> Frac_verts_host(DSIZE);
        thrust::device_vector<cuDFNsys::Fracture<_DataType_>> Frac_verts_device(DSIZE);
        cuDFNsys::Fracture<_DataType_> *Frac_verts_device_ptr;
        Frac_verts_device_ptr = thrust::raw_pointer_cast(Frac_verts_device.data());

        cuDFNsys::Warmup<<<DSIZE / 256 + 1, 256 /*  1, 2*/>>>();
        cudaDeviceSynchronize();
        time_t t;
        time(&t);

        cout << "generating fractures" << endl;

        cuDFNsys::FracturesFour<_DataType_><<<DSIZE / 256 + 1, 256 /*  1, 2*/>>>(Frac_verts_device_ptr,
                                                                                 (unsigned long)t,
                                                                                 DSIZE, L);
        cudaDeviceSynchronize();

        Frac_verts_host = Frac_verts_device;
        cout << "identifying intersections with complete fractures" << endl;
        std::map<pair<size_t, size_t>, pair<cuDFNsys::Vector3<_DataType_>, cuDFNsys::Vector3<_DataType_>>> Intersection_map;
        cuDFNsys::IdentifyIntersection<_DataType_> identifyInters{Frac_verts_host.size(),
                                                                  Frac_verts_device_ptr,
                                                                  false,
                                                                  Intersection_map};
        cout << "identifying cluster with complete fractures" << endl;
        std::vector<std::vector<size_t>> ListClusters;
        std::vector<size_t> Percolation_cluster;
        cuDFNsys::Graph<_DataType_> G{(size_t)DSIZE, Intersection_map};
        G.UseDFS(ListClusters);
        cuDFNsys::IdentifyPercolationCluster<_DataType_> IdentiClu{ListClusters,
                                                                   Frac_verts_host, perco_dir,
                                                                   Percolation_cluster};
        cout << "DFN I finished" << endl;
        cuDFNsys::MatlabPlotDFN<_DataType_> As{"DFN_I.h5", "DFN_I.m",
                                               Frac_verts_host, Intersection_map, ListClusters,
                                               Percolation_cluster, false, true, true, true,
                                               L, perco_dir, true, "DFN_I"};

        cuDFNsys::OutputObjectData<_DataType_> lk;
        // lk.OutputFractures("Fractures.h5", Frac_verts_host, L);

        //
        Intersection_map.clear();
        ListClusters.clear();
        Percolation_cluster.clear();
        cout << "identifying intersections with truncated fractures" << endl;
        cuDFNsys::IdentifyIntersection<_DataType_> identifyInters2{Frac_verts_host.size(),
                                                                   Frac_verts_device_ptr, true,
                                                                   Intersection_map};
        cout << "identifying cluster with truncated fractures" << endl;
        cuDFNsys::Graph<_DataType_> G2{(size_t)DSIZE, Intersection_map};
        G2.UseDFS(ListClusters);
        cuDFNsys::IdentifyPercolationCluster<_DataType_> IdentiClu2{ListClusters,
                                                                    Frac_verts_host, perco_dir,
                                                                    Percolation_cluster};
        cout << "DFN II finished" << endl;
        cuDFNsys::MatlabPlotDFN<_DataType_> As2{"DFN_II.h5", "DFN_II.m",
                                                Frac_verts_host, Intersection_map, ListClusters,
                                                Percolation_cluster, true, true, true, true,
                                                L, perco_dir, true, "DFN_II"};
        Frac_verts_device.clear();
        Frac_verts_device.shrink_to_fit();

        //-----------
        if (Percolation_cluster.size() > 0)
        {

            double istart_1 = cuDFNsys::CPUSecond();
            std::vector<size_t> Fracs_percol;
            cuDFNsys::GetAllPercolatingFractures GetPer{Percolation_cluster,
                                                        ListClusters,
                                                        Fracs_percol};

            std::vector<pair<int, int>> IntersectionPair_percol;

            cuDFNsys::RemoveDeadEndFrac<_DataType_> RDEF{Fracs_percol,
                                                         IntersectionPair_percol,
                                                         (size_t)perco_dir,
                                                         Frac_verts_host,
                                                         Intersection_map};
            cout << "meshing ..." << endl;

            cuDFNsys::Mesh<_DataType_> mesh{Frac_verts_host, IntersectionPair_percol,
                                            &Fracs_percol, minGrid, maxGrid, perco_dir, L};
            //lk.OutputMesh("mesh.h5", mesh, Fracs_percol);
            int i = 0;
            mesh.MatlabPlot("DFN_mesh_" + to_string(i + 1) + ".h5",
                            "N",
                            Frac_verts_host, L, true, true, true, "DFN_mesh_" + to_string(i + 1));

            cout << "MHFEM ing ..." << endl;

            cuDFNsys::MHFEM<_DataType_> fem{mesh, Frac_verts_host, 100, 20, perco_dir, L};

            cout << "Fluxes: " << fem.QIn << ", ";
            cout << fem.QOut << ", Permeability: ";
            cout << fem.Permeability << endl;
            if (fem.QError > 1 || isnan(fem.Permeability) == 1)
            {
                cout << "\e[1;32mFound large error or isnan, the error: " << fem.QError << ", the permeability: " << fem.Permeability << "\e[0m\n";
            }
            double ielaps_1 = cuDFNsys::CPUSecond() - istart_1;
            cout << "Running time of the meshing and flow simulation: ";
            cout << ielaps_1 << " sec\n";
            //---------------------
            fem.MatlabPlot("MHFEM_" + to_string(i + 1) + ".h5",
                           "MHFEM_" + to_string(i + 1) + ".m",
                           Frac_verts_host, mesh, L, true, "MHFEM_" + to_string(i + 1));
            //---------------
            lk.OutputFractures("FracturesForParticle.h5", Frac_verts_host, L);
            //return 0;
            cout << "Particle transport ing ...\n";

            int NUMparticles[LoopTime] = {100, 1000, 10000, 20000, 30000};

            int Nproc = 10;

            int MCTIME = 5;

            string filename = "CPU_check_PT_Nproc_" + std::to_string(Nproc);
            cuDFNsys::HDF5API hd5;

            string hdfilename = filename + "/" + "DFN_PT_countTime.h5";

            for (int k = 0; k < LoopTime; ++k)
            {
                cout << "step " << k + 1 << endl;
                if (k == 0)
                {
                    int us = remove(filename.c_str());

                    us = mkdir(filename.c_str(), 0777);

                    if (us == -1)
                        throw cuDFNsys::ExceptionsPause("Cannot create file " + filename);

                    hd5.NewFile(hdfilename);

                    hd5.AddDataset(hdfilename, "N",
                                   "NUMparticles", NUMparticles, make_uint2(1, LoopTime));

                    hd5.AddDataset(hdfilename, "N",
                                   "MonteCarloTimes", &MCTIME, make_uint2(1, 1));

                    hd5.AddDataset(hdfilename, "N",
                                   "LoopTime", &LoopTime, make_uint2(1, 1));
                }

                for (int j = 0; j < MCTIME; ++j)
                {
                    cuDFNsys::ParticleTransport<_DataType_> p{NUMparticles[k], // number of particle
                                                              2001,            // number of time steps
                                                              1.1813e5,        // delta T
                                                              0,               // molecular diffusion
                                                              Frac_verts_host, mesh, fem, (uint)perco_dir, -0.5f * L,
                                                              "Particle_tracking", "Flux-weighted", true, Nproc, true};

                    uint2 Dims = make_uint2(1, p.RunTimeEveryStep.size());

                    hd5.AddDataset(hdfilename, "Step_" + cuDFNsys::ToStringWithWidth(k + 1, 5),
                                   "PT_RunTime_EachStep_MC_" + cuDFNsys::ToStringWithWidth(j + 1, 5), p.RunTimeEveryStep.data(), Dims);
                }
            }
        }
        //cudaDeviceReset();
        double ielaps = cuDFNsys::CPUSecond() - istart;
        cout << "Running time of this simulation: " << ielaps << " sec\n";
    }
    catch (cuDFNsys::ExceptionsIgnore &e)
    {
        cout << e.what() << endl;
    }
    catch (cuDFNsys::ExceptionsPause &e)
    {
        cout << e.what() << endl;
    }
    catch (...)
    {
        throw;
    };
    return 0;
};