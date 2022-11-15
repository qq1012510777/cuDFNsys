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

        time_t t;
        time(&t);

        int DSIZE = 0;
        _DataType_ L = 0;
        //_DataType_ minGrid = 0;
        //_DataType_ maxGrid = 0;

        DSIZE = 256;
        L = 30;

        cuDFNsys::Warmup<<<DSIZE / 256 + 1, 256 /*  1, 2*/>>>();
        cudaDeviceSynchronize();

        int perco_dir = 2;

        cout << "preparing" << endl;

        thrust::host_vector<cuDFNsys::Fracture<_DataType_>> Frac_verts_host;
        thrust::device_vector<cuDFNsys::Fracture<_DataType_>> Frac_verts_device;

        cuDFNsys::InputObjectData<_DataType_> lk;
        lk.InputFractures("Fractures.h5", Frac_verts_host, L);
        DSIZE = Frac_verts_host.size();

        Frac_verts_device = Frac_verts_host;
        cuDFNsys::Fracture<_DataType_> *Frac_verts_device_ptr;
        Frac_verts_device_ptr = thrust::raw_pointer_cast(Frac_verts_device.data());

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
        cuDFNsys::MatlabPlotDFN<_DataType_> As{"DFN_I.h5", "N",
                                               Frac_verts_host, Intersection_map, ListClusters,
                                               Percolation_cluster, false, true, true, true,
                                               L, perco_dir, true, "DFN_I"};

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
        cuDFNsys::MatlabPlotDFN<_DataType_> As2{"DFN_II.h5", "N",
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

            cuDFNsys::OutputObjectData<_DataType_> lkew;
            lkew.OutputFractures("FracturesII.h5", Frac_verts_host, L);

            cuDFNsys::Mesh<_DataType_> mesh;

            lk.InputMesh("mesh.h5", mesh, &Fracs_percol);

            int i = 0;
            mesh.MatlabPlot("DFN_mesh_" + to_string(i + 1) + ".h5",
                            "N",
                            Frac_verts_host, L, true, true, true, "DFN_mesh_" + to_string(i + 1));

            cout << "MHFEM ing ..." << endl;

            _DataType_ P_in = 100, P_out = 20;

            if (argv[5] != NULL && argv[6] != NULL)
                P_in = atof(argv[5]), P_out = atof(argv[6]);

            cuDFNsys::MHFEM<_DataType_> fem{mesh, Frac_verts_host, P_in, P_out, perco_dir, L};
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
                           "N",
                           Frac_verts_host, mesh, L, true, "MHFEM_" + to_string(i + 1));
            //---------------
            cout << "Particle transport ing ...\n";

            cuDFNsys::OutputObjectData<_DataType_> lk;
            lk.OutputFractures("FracturesForParticle.h5", Frac_verts_host, L);

            cuDFNsys::ParticleTransport<_DataType_> p{atoi(argv[1]), // number of time step
                                                      Frac_verts_host, mesh, fem, (uint)perco_dir, -0.5f * L,
                                                      atoi(argv[2]), // num of particle
                                                      atof(argv[3]), // delta_T_ii
                                                      atof(argv[4]),
                                                      "Particle_tracking",
                                                      "Flux-weighted"};
            p.MatlabPlot("MHFEM_" + to_string(i + 1) + ".h5", "ParticlesDFNMatlab.m", mesh, fem, L, true, "ParticlesDFN");
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