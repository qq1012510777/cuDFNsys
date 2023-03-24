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
// NAME:        DispersionAtOneDensityValue.cu
// DESCRIPTION: Dispersion in a DFN with a specific percolation parameter value
// AUTHOR:      Tingchang YIN
// DATE:        24/03/2023
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
        int dev = 0;
        GPUErrCheck(cudaSetDevice(dev));
        cuDFNsys::Warmup<<<256 / 256 + 1, 256 /*  1, 2*/>>>();
        cudaDeviceSynchronize();

        time_t t;
        time(&t);

        srand((unsigned int)time(0));

        Eigen::MatrixXd Ter = Eigen::MatrixXd::Random(1, 1);

        int DSIZE = atoi(argv[1]);
        _DataType_ L = atof(argv[2]);
        _DataType_ kappa_ = atof(argv[3]),
                   beta_ = atof(argv[4]),
                   gamma_ = atof(argv[5]);
        int size_frac_mode = atoi(argv[6]); // mode of fracture size distributions
        cuDFNsys::Vector4<_DataType_> ParaSizeDistri =
            cuDFNsys::MakeVector4((_DataType_)atof(argv[7]),
                                  (_DataType_)atof(argv[8]),
                                  (_DataType_)atof(argv[9]),
                                  (_DataType_)atof(argv[10]));
        double3 DomainDimensionRatio = make_double3(1, 1, atof(argv[11]));
        int IfRemoveDeadEnd = atoi(argv[12]);
        _DataType_ minGridSize = atof(argv[13]);
        _DataType_ maxGridSize = atof(argv[14]);
        _DataType_ P_in = atof(argv[15]), P_out = atof(argv[16]);

        int NumTimeSteps_Dispersion = atoi(argv[17]);
        int NumParticlesRandomWalk = atoi(argv[18]);
        _DataType_ DeltaT = 0;
        _DataType_ DiffusionLocal = atof(argv[19]);
        _DataType_ ControlPlaneSpacing = atof(argv[20]);

        cout << "Number of fractures: " << DSIZE << endl;
        cout << "L: " << L << endl;
        cout << "Kappa: " << kappa_ << endl;
        cout << "Beta: " << beta_ << endl;
        cout << "Gamma: " << gamma_ << endl;
        cout << "Mode of fracture size distributions: " << size_frac_mode << endl;
        cout << "Parameters of the size distribution: " << ParaSizeDistri.x << ", " << ParaSizeDistri.y << ", " << ParaSizeDistri.z << ", " << ParaSizeDistri.w << endl;
        cout << "Domain's dimension ratio: " << DomainDimensionRatio.x << ", " << DomainDimensionRatio.y << ", " << DomainDimensionRatio.z << endl;
        cout << "If remove the dead ends: " << (IfRemoveDeadEnd == 0 ? "false" : "true") << endl;
        cout << "Min grid size: " << minGridSize << endl;
        cout << "Max grid size: " << maxGridSize << endl;
        cout << "Hydraulic head at the inlet and outlet: " << P_in << ", " << P_out << endl;
        cout << "Number of time steps for random walks: " << NumTimeSteps_Dispersion << endl;
        cout << "Number of particles: " << NumParticlesRandomWalk << endl;
        //cout << "Delta T: " << DeltaT << endl;
        cout << "Molecular diffusion: " << DiffusionLocal << endl;
        cout << "The spacing of control planes: " << ControlPlaneSpacing << endl;

        int perco_dir = 2;

        string FractureFileName = "Fractures.h5";

        std::ifstream fileer(FractureFileName);
        bool pwqsc = fileer.good();

        if (!pwqsc) // no DFN is existing
        {
            thrust::host_vector<cuDFNsys::Fracture<_DataType_>> Frac_verts_host(DSIZE);
            thrust::device_vector<cuDFNsys::Fracture<_DataType_>> Frac_verts_device(DSIZE);
            cuDFNsys::Fracture<_DataType_> *Frac_verts_device_ptr;
            Frac_verts_device_ptr = thrust::raw_pointer_cast(Frac_verts_device.data());

            cuDFNsys::Fractures<_DataType_><<<DSIZE / 256 + 1, 256>>>(Frac_verts_device_ptr,
                                                                      (unsigned long)t + (unsigned long)ceil(abs(Ter(0, 0)) * ((unsigned long)t * 1.0)),
                                                                      DSIZE, L,
                                                                      0,
                                                                      ParaSizeDistri,
                                                                      kappa_, // kappa
                                                                      beta_,  // beta
                                                                      gamma_, // gamma
                                                                      DomainDimensionRatio);
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
                // double istart_1 = cuDFNsys::CPUSecond();
                std::vector<size_t> Fracs_percol;

                cuDFNsys::GetAllPercolatingFractures GetPer{Percolation_cluster,
                                                            ListClusters,
                                                            Fracs_percol};
                std::vector<pair<int, int>> IntersectionPair_percol;
                int NUMprior = Fracs_percol.size();

                bool ifRemoveDeadEnds = (IfRemoveDeadEnd == 0 ? false : true);
                cout << "ifRemoveDeadEnds: " << ifRemoveDeadEnds << endl;
                cuDFNsys::RemoveDeadEndFrac<_DataType_> RDEF{Fracs_percol,
                                                             IntersectionPair_percol,
                                                             (size_t)perco_dir,
                                                             Frac_verts_host,
                                                             Intersection_map,
                                                             ifRemoveDeadEnds};

                if (ifRemoveDeadEnds)
                    cout << "remove " << NUMprior - Frac_verts_host.size() << " fractures\n";

                cout << "meshing ..." << endl;

                cuDFNsys::Mesh<_DataType_> mesh{Frac_verts_host, IntersectionPair_percol,
                                                &Fracs_percol, minGridSize, maxGridSize, perco_dir, L,
                                                DomainDimensionRatio};
                lk.OutputMesh("mesh.h5", mesh, Fracs_percol);

                mesh.MatlabPlot("DFN_mesh_.h5",
                                "DFN_mesh_.m",
                                Frac_verts_host,
                                L, true, true, true,
                                "DFN_mesh_", DomainDimensionRatio);

                cout << "MHFEM ing ..." << endl;

                cuDFNsys::MHFEM<_DataType_> fem{mesh, Frac_verts_host, P_in, P_out,
                                                perco_dir, L, DomainDimensionRatio};
                lk.OutputMHFEM("mhfem.h5", fem);

                cout << "Fluxes: " << fem.QIn << ", ";
                cout << fem.QOut << ", Permeability: ";
                cout << fem.Permeability << endl;
                cout << "Error between the inlet and outlet fluxes: " << abs(fem.QIn - fem.QOut) / ((fem.QOut + fem.QIn) * 0.5) * 100.0 << "%\n";
                if (fem.QError > 1 || isnan(fem.Permeability) == 1)
                    throw cuDFNsys::ExceptionsIgnore("Found large error or isnan, the error: " + std::to_string(fem.QError) + ", the permeability: " + std::to_string(fem.Permeability) + "\n");

                //---------------------
                double meanV = fem.MatlabPlot("MHFEM_.h5",
                                              "MHFEM_.m",
                                              Frac_verts_host, mesh, L, true, "MHFEM_", DomainDimensionRatio);
                cout << "The mean velocity of all elements is " << meanV << endl;

                double ik = 0;
                for (int u = 1;; u++)
                    if (pow(10, u) * meanV >= 1.0)
                    {
                        ik = u;
                        break;
                    }

                DeltaT = 0.5 * pow(10, ik - 1);

                cout << "\nThe delta T is set to be " << ("\033[1;33m") << DeltaT << ("\033[0m") << "\n\n";
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

                cout << "Particle transport ing ......\n";

                cuDFNsys::ParticleTransport<_DataType_> p{NumTimeSteps_Dispersion, // number of time step
                                                          Frac_verts_host,
                                                          mesh,
                                                          fem,
                                                          (uint)perco_dir,
                                                          -0.5f * L * (&DomainDimensionRatio.x)[perco_dir],
                                                          NumParticlesRandomWalk, // num of particle
                                                          DeltaT,                 // delta_T_ii
                                                          DiffusionLocal,
                                                          "Particle_tracking",
                                                          "Flux-weighted",
                                                          "FPTCurve",
                                                          false, 1, false,
                                                          ControlPlaneSpacing};

                p.MatlabPlot("MHFEM_.h5", "ParticlesDFNMatlab.m", mesh, fem, L, DomainDimensionRatio, true, "ParticlesDFN");
            }
        }
        else // if the file exists
        {
            cout << "preparing" << endl;

            thrust::host_vector<cuDFNsys::Fracture<_DataType_>> Frac_verts_host;
            thrust::device_vector<cuDFNsys::Fracture<_DataType_>> Frac_verts_device;

            cuDFNsys::InputObjectData<_DataType_> lk;
            lk.InputFractures("Fractures.h5", Frac_verts_host, L, DomainDimensionRatio);

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
            cuDFNsys::MatlabPlotDFN<_DataType_> As{"DFN_I.h5", "DFN_I.m",
                                                   Frac_verts_host, Intersection_map, ListClusters,
                                                   Percolation_cluster, false, true, true, true,
                                                   L, perco_dir, true, "DFN_I", DomainDimensionRatio};

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

                bool ifRemoveDeadends = (IfRemoveDeadEnd == 0 ? false : true);

                cuDFNsys::RemoveDeadEndFrac<_DataType_> RDEF{Fracs_percol,
                                                             IntersectionPair_percol,
                                                             (size_t)perco_dir,
                                                             Frac_verts_host,
                                                             Intersection_map, ifRemoveDeadends};
                cout << "meshing ..." << endl;

                cuDFNsys::OutputObjectData<_DataType_> lkew;
                lkew.OutputFractures("FracturesII.h5", Frac_verts_host, L, DomainDimensionRatio);

                cuDFNsys::Mesh<_DataType_> mesh;
                try
                {
                    lk.InputMesh("mesh.h5", mesh, &Fracs_percol);
                }
                catch (...)
                {
                    cout << "mesh ing ...\n";
                    cuDFNsys::Mesh<_DataType_> mesh2{Frac_verts_host, IntersectionPair_percol,
                                                     &Fracs_percol, minGridSize,
                                                     maxGridSize, perco_dir, L, DomainDimensionRatio};
                    lkew.OutputMesh("mesh.h5", mesh2, Fracs_percol);
                    lk.InputMesh("mesh.h5", mesh, &Fracs_percol);
                }

                mesh.MatlabPlot("DFN_mesh_.h5",
                                "DFN_mesh_.m",
                                Frac_verts_host, L, true, true, true, "DFN_mesh_", DomainDimensionRatio);

                cout << "MHFEM ing ..." << endl;

                cuDFNsys::MHFEM<_DataType_> fem;
                try
                {
                    cout << "Loading mhfem ...\n";
                    lk.InputMHFEM("mhfem.h5", fem);
                }
                catch (...)
                {
                    cuDFNsys::MHFEM<_DataType_> fem2{mesh, Frac_verts_host, P_in, P_out, perco_dir, L, DomainDimensionRatio};

                    lkew.OutputMHFEM("mhfem.h5", fem2);
                    fem = fem2;
                };

                cout << "Fluxes: " << fem.QIn << ", ";
                cout << fem.QOut << ", Permeability: ";
                cout << fem.Permeability << endl;
                if (fem.QError > 1 || isnan(fem.Permeability) == 1)
                    throw cuDFNsys::ExceptionsIgnore("Found large error or isnan, the error: " + std::to_string(fem.QError) + ", the permeability: " + std::to_string(fem.Permeability) + "\n");

                //---------------------
                double meanV = fem.MatlabPlot("MHFEM_.h5",
                                              "MHFEM_.m",
                                              Frac_verts_host, mesh, L, true, "MHFEM_", DomainDimensionRatio);
                cout << "The mean velocity of all elements is " << meanV << endl;

                double ik = 0;
                for (int u = 1;; u++)
                    if (pow(10, u) * meanV >= 1.0)
                    {
                        ik = u;
                        break;
                    }

                DeltaT = 0.5 * pow(10, ik - 1);

                cout << "\nThe delta T is set to be " << ("\033[1;33m") << DeltaT << ("\033[0m") << "\n\n";
                //---------------
                // return 0;

                string Filename_FracturesForParticle = "FracturesForParticle.h5";

                std::ifstream file(Filename_FracturesForParticle);
                bool pwqs = file.good();

                if (!pwqs)
                {
                    cout << "Writting " << Filename_FracturesForParticle << endl;
                    cuDFNsys::OutputObjectData<_DataType_> lk;
                    lk.OutputFractures(Filename_FracturesForParticle, Frac_verts_host, L);
                }
                cout << "Particle transport ing ......\n";

                cuDFNsys::ParticleTransport<_DataType_> p{NumTimeSteps_Dispersion, // number of time step
                                                          Frac_verts_host, mesh, fem, (uint)perco_dir,
                                                          -0.5f * L * (&DomainDimensionRatio.x)[perco_dir],
                                                          NumParticlesRandomWalk, // num of particle
                                                          DeltaT,                 // delta_T_ii
                                                          DiffusionLocal,
                                                          "Particle_tracking",
                                                          "Flux-weighted",
                                                          "FPTCurve",
                                                          false, 1, false, ControlPlaneSpacing};
                p.MatlabPlot("MHFEM_.h5", "ParticlesDFNMatlab.m", mesh, fem, L, DomainDimensionRatio, true, "ParticlesDFN");
            }
            //cudaDeviceReset();
        }
    }
    catch (cuDFNsys::ExceptionsIgnore &e)
    {
        cout << e.what() << endl;
        exit(0);
    }
    catch (cuDFNsys::ExceptionsPause &e)
    {
        cout << e.what() << endl;
        exit(0);
    }
    catch (...)
    {
        throw;
    }
    return 0;
};
