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
// NAME:        IdenticalDFNDifferentPe.cu
// DESCRIPTION: 1
// AUTHOR:      Tingchang YIN
// DATE:        12/06/2023
// ====================================================

#include "cuDFNsys.cuh"
#include <H5Exception.h>
#include <fstream>
#include <iostream>
#include <limits.h>
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
        _DataType_ Factor_mean_time_in_grid = atof(argv[1]);
        _DataType_ LengthScale_Over_Pe = atof(argv[2]) / atof(argv[3]);
        int NumTimeSteps_Dispersion = atoi(argv[4]);
        int NumParticlesRandomWalk = atoi(argv[5]);
        string injectionMode = argv[6]; // Flux-weighted or Resident
        string OutputMode = argv[7];    // OutputAll or FPTCurve
        bool IfInjectAt_Center = (atoi(argv[8]) == 0 ? false : true);
        _DataType_ InjectionPlane = atof(argv[9]); 

        int perco_dir = 2;

        std::ifstream fileer1("Fractures.h5");
        std::ifstream fileer2("mesh.h5");
        std::ifstream fileer3("mhfem.h5");

        if (fileer1.good() && fileer2.good() && fileer3.good())
        {
            // do nothing
        }
        else
        {
            cout << "No DFN, mesh, and mhfem files\n";
            exit(0);
        }

        int dev = 0;
        GPUErrCheck(cudaSetDevice(dev));
        cuDFNsys::Warmup<<<256 / 256 + 1, 256 /*  1, 2*/>>>();
        cudaDeviceSynchronize();

        time_t t;
        time(&t);

        srand((unsigned int)time(0));

        Eigen::MatrixXd Ter = Eigen::MatrixXd::Random(1, 1);

        thrust::host_vector<cuDFNsys::Fracture<_DataType_>> Frac_verts_host;
        thrust::device_vector<cuDFNsys::Fracture<_DataType_>> Frac_verts_device;
        cuDFNsys::Fracture<_DataType_> *Frac_verts_device_ptr;

        cuDFNsys::InputObjectData<_DataType_> lk;

        Frac_verts_host.resize(0);

        double3 DomainDimensionRatio;
        _DataType_ L;

        lk.InputFractures("Fractures.h5", Frac_verts_host, L, DomainDimensionRatio);
        Frac_verts_host.shrink_to_fit();
        Frac_verts_device = Frac_verts_host;
        Frac_verts_device_ptr = thrust::raw_pointer_cast(Frac_verts_device.data());
        uint DSIZE = Frac_verts_host.size();

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

        // std::vector<double> Data1(14);
        // cuDFNsys::GetStatistics<double>(Frac_verts_host,
        //                                 Intersection_map,
        //                                 ListClusters,
        //                                 Percolation_cluster, L, Data1[0], Data1[1], Data1[2], Data1[3],
        //                                 Data1[4], Data1[5], Data1[6], Data1[7], Data1[8], Data1[9],
        //                                 Data1[10], Data1[11], Data1[12], Data1[13]);
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
                                                L, perco_dir, true, "DFN_II", DomainDimensionRatio};
        //return 0;
        Frac_verts_device.clear();
        Frac_verts_device.shrink_to_fit();

        cuDFNsys::OutputObjectData<_DataType_> lk_out;

        if (Percolation_cluster.size() > 0)
        {
            //lk_out.OutputFractures("Fractures.h5", Frac_verts_host, L, DomainDimensionRatio);

            std::vector<size_t> Fracs_percol;

            cuDFNsys::GetAllPercolatingFractures GetPer{Percolation_cluster,
                                                        ListClusters,
                                                        Fracs_percol};
            std::vector<pair<int, int>> IntersectionPair_percol;
            int NUMprior = Fracs_percol.size();

            bool ifRemoveDeadEnds = false;
            cout << "ifRemoveDeadEnds: " << ifRemoveDeadEnds << endl;
            cuDFNsys::RemoveDeadEndFrac<_DataType_> RDEF{Fracs_percol,
                                                         IntersectionPair_percol,
                                                         (size_t)perco_dir,
                                                         Frac_verts_host,
                                                         Intersection_map,
                                                         ifRemoveDeadEnds};

            cout << "meshing ..." << endl;

            cuDFNsys::Mesh<_DataType_> mesh;

            lk_out.OutputFractures("FracturesII.h5", Frac_verts_host, L, DomainDimensionRatio);

            lk.InputMesh("mesh.h5", mesh, &Fracs_percol);

            double mean_grid_area = mesh.MatlabPlot("DFN_mesh_.h5",
                                                    "DFN_mesh_.m",
                                                    Frac_verts_host,
                                                    L, true, true, true,
                                                    "DFN_mesh_", DomainDimensionRatio);

            cout << "The mean area of all elements is " << mean_grid_area << endl;

            cout << "MHFEM ing ..." << endl;

            cuDFNsys::MHFEM<_DataType_> fem;

            lk.InputMHFEM("mhfem.h5", fem);

            cout << "Fluxes: " << fem.QIn << ", ";
            cout << fem.QOut << ", Permeability: ";
            cout << fem.Permeability << endl;
            if (fem.QError > 1 || isnan(fem.Permeability) == 1)
            {
                throw cuDFNsys::ExceptionsIgnore("Found large error or isnan, the error: " + std::to_string(fem.QError) + ", the permeability: " + std::to_string(fem.Permeability) + "\n");
            }
            double2 TGH = fem.MatlabPlot("MHFEM_.h5",
                                         "MHFEM_.m",
                                         Frac_verts_host, mesh, L, true, "MHFEM_", DomainDimensionRatio);

            double meanV = TGH.x;
            double maxV = TGH.y;

            cout << "The maximum velocity of all elements is " << maxV << endl;
            cout << "The mean velocity of all elements is " << meanV << endl;

            double DeltaT;
            double DiffusionLocal;

            double meanTime = pow(mean_grid_area, 0.5) / maxV;

            DeltaT = meanTime / Factor_mean_time_in_grid;
            cout << "\nThe delta T is set to be " << ("\033[1;33m") << DeltaT << ("\033[0m") << "\n\n";

            DiffusionLocal = LengthScale_Over_Pe * meanV;
            cout << "\nThe DiffusionLocal is set to be " << ("\033[1;33m") << DiffusionLocal << ("\033[0m") << "\n\n";
            //-----------------

            //string FractureFileName_r = "ParticlePositionResult/DispersionInfo.h5";

            string Filename_FracturesForParticle = "FracturesForParticle.h5";

            std::ifstream file(Filename_FracturesForParticle);
            bool pwqs = file.good();

            if (!pwqs)
            {
                cout << "Writting " << Filename_FracturesForParticle << endl;
                lk_out.OutputFractures(Filename_FracturesForParticle, Frac_verts_host, L);
            }

            cout << "Particle transport ing ......\n";

            cuDFNsys::ParticleTransport<_DataType_> p{NumTimeSteps_Dispersion, // number of time step
                                                      Frac_verts_host, mesh, fem, (uint)perco_dir,
                                                      -0.5f * L * (&DomainDimensionRatio.x)[perco_dir],
                                                      NumParticlesRandomWalk, // num of particle
                                                      DeltaT,                 // delta_T_ii
                                                      DiffusionLocal,
                                                      "Particle_tracking",
                                                      injectionMode, // Flux-weighted or Resident
                                                      OutputMode,
                                                      false, 1, false, 10000, true, IfInjectAt_Center, InjectionPlane};
            p.MatlabPlot("MHFEM_.h5", "ParticlesDFNMatlab.m", mesh, fem, L, DomainDimensionRatio, true, "ParticlesDFN");
        }
    }
    catch (cuDFNsys::ExceptionsIgnore &e)
    {
        cout << "cuDFNsys::ExceptionsIgnore\n";
        cout << e.what() << endl;
    }
    catch (cuDFNsys::ExceptionsPause &e)
    {
        cout << "cuDFNsys::ExceptionsPause\n";
        cout << e.what() << endl;
    }
    catch (...)
    {
        throw;
    }
    return 0;
};
