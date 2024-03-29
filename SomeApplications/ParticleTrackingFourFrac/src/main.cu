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
        /// uint iyu[1000] = {1};
        /// uint2 dim_s = make_uint2(1, 1000);
        /// cuDFNsys::HDF5API hg0;
        /// vector<string> datasetname = {"kl", "lio"};
        /// vector<uint *> sd = {iyu, iyu};
        /// vector<uint2> ee = {dim_s, dim_s};
        /// hg0.NewFile("Test1.h5");
        /// hg0.AddDatasetsWithOneGroup("Test1.h5", "N", datasetname, sd, ee);
        /// vector<uint> sdd = hg0.ReadDataset<uint>("Test1.h5", "N", "lio");
        /// cout << sdd[0] << ", " << sdd[1] << endl;
        /// return 0;

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
            // return 0;
            cout << "Particle transport ing ...\n";

            int NumTimeStep = atoi(argv[1]);
            int NumParticles = atoi(argv[2]);
            _DataType_ DeltaT = 1.1813e5;
            _DataType_ DiffusionMole = 0.0;

            double3 DomainDimensionRatio = make_double3(1, 1, 1);

            cuDFNsys::ParticleTransport<_DataType_> p{NumTimeStep,                                     // number of time steps
                                                      Frac_verts_host,                                 // fractures
                                                      mesh,                                            // mesh
                                                      fem,                                             // fem
                                                      (uint)perco_dir,                                 // percolation direction
                                                      -0.5 * (&DomainDimensionRatio.x)[perco_dir] * L, // the target plane, z = -0.5 * L; it could be changed if the percolation direction is not along the z axis
                                                      NumParticles,                                    // number of particles
                                                      DeltaT,                                          // delta t
                                                      DiffusionMole,                                   // molecular diffusion
                                                      "Particle_tracking",                             // use particle tracking algorithm
                                                      "Flux-weighted",                                 // injection mode
                                                      "OutputAll",                                     // output all information
                                                      false,                                           // if use CPU?
                                                      0,                                               // number of CPU processors
                                                      false,                                           // if record run time
                                                      50};                                             // the spacing of the linear distance of control planes

            p.MatlabPlot("MHFEM_" + to_string(i + 1) + ".h5", "ParticleTracking.m", mesh, fem, L, DomainDimensionRatio, true, "particlePy");
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