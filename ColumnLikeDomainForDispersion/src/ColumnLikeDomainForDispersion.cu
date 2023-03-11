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
// NAME:        ColumnLikeDomainForDispersion.cu
// DESCRIPTION: column-like domain
// AUTHOR:      Tingchang YIN
// DATE:        10/03/2023
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
        cuDFNsys::Warmup<<<256 / 256 + 1, 256 /*  1, 2*/>>>();
        cudaDeviceSynchronize();

        time_t t;
        time(&t);

        srand((unsigned int)time(0));
        Eigen::MatrixXd Ter = Eigen::MatrixXd::Random(1, 1);
        int DSIZE = 800;
        _DataType_ L = 100;
        cuDFNsys::Vector4<_DataType_> ParaSizeDistri =
            cuDFNsys::MakeVector4((_DataType_)1.5,
                                  (_DataType_)1,
                                  (_DataType_)100,
                                  (_DataType_)0);
        double3 DomainDimensionRatio = make_double3(1, 1, 4);
        int perco_dir = 2;

        thrust::host_vector<cuDFNsys::Fracture<_DataType_>> Frac_verts_host(DSIZE);
        thrust::device_vector<cuDFNsys::Fracture<_DataType_>> Frac_verts_device(DSIZE);
        cuDFNsys::Fracture<_DataType_> *Frac_verts_device_ptr;
        Frac_verts_device_ptr = thrust::raw_pointer_cast(Frac_verts_device.data());

        cuDFNsys::Fractures<_DataType_><<<DSIZE / 256 + 1, 256>>>(Frac_verts_device_ptr,
                                                                  (unsigned long)t + (unsigned long)ceil(abs(Ter(0, 0)) * ((unsigned long)t * 1.0)),
                                                                  DSIZE, L,
                                                                  0,
                                                                  ParaSizeDistri,
                                                                  0,                             // kappa
                                                                  0.1,                           // beta
                                                                  5.2e-4, DomainDimensionRatio); // gamma
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
                                               L, perco_dir, true, "DFN_I", DomainDimensionRatio};
        cuDFNsys::OutputObjectData<_DataType_> lk;
        lk.OutputFractures("Fractures.h5", Frac_verts_host, L, DomainDimensionRatio);

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
                                                L, perco_dir, true, "DFN_II", DomainDimensionRatio};
        //return 0;
        Frac_verts_device.clear();
        Frac_verts_device.shrink_to_fit();

        if (Percolation_cluster.size() > 0)
        {

            double istart_1 = cuDFNsys::CPUSecond();
            std::vector<size_t> Fracs_percol;

            cuDFNsys::GetAllPercolatingFractures GetPer{Percolation_cluster,
                                                        ListClusters,
                                                        Fracs_percol};
            std::vector<pair<int, int>> IntersectionPair_percol;
            int NUMprior = Fracs_percol.size();

            bool ifRemoveDeadEnds = (atoi(argv[12]) == 0 ? false : true);

            cuDFNsys::RemoveDeadEndFrac<_DataType_> RDEF{Fracs_percol,
                                                         IntersectionPair_percol,
                                                         (size_t)perco_dir,
                                                         Frac_verts_host,
                                                         Intersection_map, ifRemoveDeadEnds};

            if (ifRemoveDeadEnds)
                cout << "remove " << NUMprior - Frac_verts_host.size() << " fractures\n";

            cout << "meshing ..." << endl;

            cuDFNsys::Mesh<_DataType_> mesh{Frac_verts_host, IntersectionPair_percol,
                                            &Fracs_percol, 1, 10, perco_dir, L,
                                            DomainDimensionRatio};
            lk.OutputMesh("mesh.h5", mesh, Fracs_percol);

            mesh.MatlabPlot("DFN_mesh_.h5",
                            "DFN_mesh_.m",
                            Frac_verts_host,
                            L, true, true, true,
                            "DFN_mesh_", DomainDimensionRatio);

            cout << "MHFEM ing ..." << endl;

            cuDFNsys::MHFEM<_DataType_> fem{mesh, Frac_verts_host, 100, 20,
                                            perco_dir, L, DomainDimensionRatio};
            lk.OutputMHFEM("mhfem.h5", fem);

            cout << "Fluxes: " << fem.QIn << ", ";
            cout << fem.QOut << ", Permeability: ";
            cout << fem.Permeability << endl;
            cout << "Error between the inlet and outlet fluxes: " << abs(fem.QIn - fem.QOut) / ((fem.QOut + fem.QIn) * 0.5) * 100.0 << "%\n";
            if (fem.QError > 1 || isnan(fem.Permeability) == 1)
            {
                cout << "\e[1;32mFound large error or isnan, the error: " << fem.QError << ", the permeability: " << fem.Permeability << "\e[0m\n";
            }

            //---------------------
            fem.MatlabPlot("MHFEM_.h5",
                           "MHFEM_.m",
                           Frac_verts_host, mesh, L, true, "MHFEM_", DomainDimensionRatio);
            //---------------

            string Filename_FracturesForParticle = "FracturesForParticle.h5";

            std::ifstream file(Filename_FracturesForParticle);
            bool pwqs = file.good();

            if (!pwqs)
            {
                cout << "Writting " << Filename_FracturesForParticle << endl;
                cuDFNsys::OutputObjectData<_DataType_> lk;
                lk.OutputFractures(Filename_FracturesForParticle, Frac_verts_host, L, DomainDimensionRatio);
            }
            return 0;
            
            cout << "Particle transport ing ...\n";

            //double *FG = &DomainDimensionRatio.x;

            cuDFNsys::ParticleTransport<_DataType_> p{atoi(argv[1]),             // number of particle
                                                      atoi(argv[2]),             // number of time steps
                                                      (_DataType_)atof(argv[3]), // delta T
                                                      (_DataType_)atof(argv[4]), // molecular diffusion
                                                      Frac_verts_host, mesh, fem, (uint)perco_dir,
                                                      -0.5f * L * (&DomainDimensionRatio.x)[perco_dir],
                                                      "Particle_tracking",
                                                      "Flux-weighted"};

            p.MatlabPlot("MHFEM_.h5", "particle.m", mesh, fem, L, DomainDimensionRatio);
        }
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