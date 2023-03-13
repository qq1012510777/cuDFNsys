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
// NAME:        benchmark cases
// DESCRIPTION: Call cuDFNsys functions to do simulation.
// AUTHOR:      Tingchang YIN
// DATE:        20/04/2022
// ====================================================

#include "cuDFNsys.cuh"
#include <unistd.h>

#ifdef USE_DOUBLES
typedef double _DataType_;
#else
typedef float _DataType_;
#endif

int main(int argc, char *argv[])
{
    try
    {
        double istart = cuDFNsys::CPUSecond();

        int dev = 0;
        GPUErrCheck(cudaSetDevice(dev));

        int DSIZE = 0;
        _DataType_ L = 0;
        _DataType_ minGrid = 0;
        _DataType_ maxGrid = 0;

        L = 30;
        DSIZE = 1;
        minGrid = 1;
        maxGrid = 5;

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
        cuDFNsys::FracturesCrossedVertical<_DataType_><<<DSIZE / 256 + 1, 256>>>(Frac_verts_device_ptr,
                                                                                 (unsigned long)t,
                                                                                 DSIZE,
                                                                                 L);

        cudaDeviceSynchronize();

        Frac_verts_host = Frac_verts_device;

        std::map<pair<size_t, size_t>, pair<cuDFNsys::Vector3<_DataType_>, cuDFNsys::Vector3<_DataType_>>> Intersection_map;
        // cout << "Frac_verts_host.size(): " << Frac_verts_host.size() << endl;
        // cout << Frac_verts_host[0].Verts3D[0].x << ", " << Frac_verts_host[0].Verts3D[0].y << ", " << Frac_verts_host[0].Verts3D[0].z << endl;
        // cout << Frac_verts_host[0].Verts3D[1].x << ", " << Frac_verts_host[0].Verts3D[1].y << ", " << Frac_verts_host[1].Verts3D[0].z << endl;
        // cout << Frac_verts_host[0].Verts3D[2].x << ", " << Frac_verts_host[0].Verts3D[2].y << ", " << Frac_verts_host[2].Verts3D[0].z << endl;
        // cout << Frac_verts_host[0].Verts3D[3].x << ", " << Frac_verts_host[0].Verts3D[3].y << ", " << Frac_verts_host[3].Verts3D[0].z << endl;

        // cout << "identifying intersections with complete fractures" << endl;
        // cuDFNsys::IdentifyIntersection<_DataType_> identifyInters{Frac_verts_host.size(),
        //                                                           Frac_verts_device_ptr,
        //                                                           false,
        //                                                           Intersection_map};
        // cout << "identifying cluster with complete fractures" << endl;
        std::vector<std::vector<size_t>> ListClusters(1);
        ListClusters[0] = std::vector<size_t>(0);

        std::vector<size_t> Percolation_cluster(1);
        Percolation_cluster[0] = 0;

        cuDFNsys::MatlabPlotDFN<_DataType_> As2{"DFN_II.h5", "DFN_II.m",
                                                Frac_verts_host, Intersection_map, ListClusters,
                                                Percolation_cluster, true, true, true, true,
                                                L, perco_dir, true, "DFN_II"};

        Frac_verts_device.clear();
        Frac_verts_device.shrink_to_fit();
        //-----------
        cuDFNsys::OutputObjectData<_DataType_> lk;
        if (Percolation_cluster.size() > 0)
        {

            double istart_1 = cuDFNsys::CPUSecond();
            std::vector<size_t> Fracs_percol;

            // cuDFNsys::GetAllPercolatingFractures GetPer{Percolation_cluster,
            //                                             ListClusters,
            //                                             Fracs_percol};
            // std::vector<pair<int, int>> IntersectionPair_percol;
            //
            // cuDFNsys::RemoveDeadEndFrac<_DataType_> RDEF{Fracs_percol,
            //                                              IntersectionPair_percol,
            //                                              (size_t)perco_dir,
            //                                              Frac_verts_host,
            //                                              Intersection_map};

            Fracs_percol = Percolation_cluster;
            std::vector<pair<int, int>> IntersectionPair_percol(0);
            cout << "meshing ..." << endl;

            cuDFNsys::Mesh<_DataType_> mesh{Frac_verts_host, IntersectionPair_percol,
                                            &Fracs_percol, minGrid, maxGrid, perco_dir, L};

            lk.OutputMesh("mesh.h5", mesh, Fracs_percol);

            //lk.OutputFractures("FracturesII.h5", Frac_verts_host, L, make_double3(1, 1, 1));

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
            int i = 0;
            mesh.MatlabPlot("DFN_mesh_" + to_string(i + 1) + ".h5",
                            "N",
                            Frac_verts_host, L, true, true, true, "DFN_mesh_" + to_string(i + 1));
            fem.MatlabPlot("MHFEM_" + to_string(i + 1) + ".h5",
                           "MHFEM_" + to_string(i + 1) + ".m",
                           Frac_verts_host, mesh, L, true, "MHFEM_" + to_string(i + 1));
            //---------------

            string Filename_FracturesForParticle = "FracturesForParticle.h5";

            std::ifstream file(Filename_FracturesForParticle);
            bool pwqs = file.good();

            if (!pwqs)
            {
                cout << "Writting " << Filename_FracturesForParticle << endl;
                cuDFNsys::OutputObjectData<_DataType_> lk;
                lk.OutputFractures(Filename_FracturesForParticle, Frac_verts_host, L);
            }

            cout << "Particle transport ing ...\n";
            // return 0;
            cuDFNsys::ParticleTransport<_DataType_> p{100000,  // number of particle
                                                      5000,    // number of time steps
                                                      5e4,     // delta T
                                                      2.22e-8, // molecular diffusion
                                                      Frac_verts_host, mesh, fem, (uint)perco_dir, -0.5f * L,
                                                      "Particle_tracking", "Flux-weighted",
                                                      false, 1, true, "FPTCurve"};
            p.MatlabPlot("MHFEM_" + to_string(i + 1) + ".h5", "particle.m", mesh, fem, L);
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