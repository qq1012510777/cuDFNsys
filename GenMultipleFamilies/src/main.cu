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
// NAME:        main.cu
// DESCRIPTION: Generate multiple families.
// AUTHOR:      Tingchang YIN
// DATE:        24/02/2023
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
        int perco_dir = 2;
        _DataType_ L = atof(argv[1]);

        uint NumFamilies = atoi(argv[2]);
        uint count_Argv_index = 2;

        thrust::host_vector<cuDFNsys::Fracture<_DataType_>> Frac_verts_host;

        for (uint i = 0; i < NumFamilies; ++i)
        {
            uint DSIZE = atoi(argv[count_Argv_index + 1]);
            count_Argv_index++;

            cuDFNsys::Vector4<_DataType_> ParaSizeDistri =
                cuDFNsys::MakeVector4((_DataType_)atof(argv[count_Argv_index + 1]),
                                      (_DataType_)atof(argv[count_Argv_index + 2]),
                                      (_DataType_)atof(argv[count_Argv_index + 3]),
                                      (_DataType_)atof(argv[count_Argv_index + 4]));
            count_Argv_index += 4;

            _DataType_ kappa_ = (_DataType_)atof(argv[count_Argv_index + 1]);
            _DataType_ beta_ = (_DataType_)atof(argv[count_Argv_index + 2]);
            _DataType_ gamma_ = (_DataType_)atof(argv[count_Argv_index + 3]);
            count_Argv_index += 3;

            cuDFNsys::Vector3<_DataType_> rotationAxis = cuDFNsys::MakeVector3(
                (_DataType_)atof(argv[count_Argv_index + 1]),
                (_DataType_)atof(argv[count_Argv_index + 2]),
                (_DataType_)atof(argv[count_Argv_index + 3]));
            count_Argv_index += 3;

            _DataType_ AngleRotation =
                (_DataType_)atof(argv[count_Argv_index + 1]); // degree
            count_Argv_index++;

            thrust::host_vector<cuDFNsys::Fracture<_DataType_>> Frac_verts_host_i(DSIZE);
            thrust::device_vector<cuDFNsys::Fracture<_DataType_>> Frac_verts_device_i(DSIZE);
            cuDFNsys::Fracture<_DataType_> *Frac_verts_device_ptr_i;
            Frac_verts_device_ptr_i = thrust::raw_pointer_cast(Frac_verts_device_i.data());

            srand((unsigned int)time(0) + i * 100);
            Eigen::MatrixXd Ter = Eigen::MatrixXd::Random(1, 1);

            cout << "generating fracture group No. " << i + 1 << ", random seed: " << (unsigned long)ceil(abs(Ter(0, 0)) * ((unsigned long)t * 1.0)) << endl;
            cout << "\tnum of fractures: " << DSIZE << endl;
            cout << "\tangleRotation: " << AngleRotation * M_PI / 180.0 << endl;
            cout << "\trotationAxis: " << rotationAxis.x << ", " << rotationAxis.y << ", " << rotationAxis.z << endl;
            cuDFNsys::Fractures<_DataType_><<<DSIZE / 256 + 1, 256>>>(Frac_verts_device_ptr_i,
                                                                      (unsigned long)t + (unsigned long)ceil(abs(Ter(0, 0)) * ((unsigned long)t * 1.0)),
                                                                      DSIZE, L,
                                                                      0, // power law
                                                                      ParaSizeDistri,
                                                                      kappa_,  // kappa
                                                                      beta_,   // beta
                                                                      gamma_); // gamma
            cudaDeviceSynchronize();
            Frac_verts_host_i = Frac_verts_device_i;

            //-------------rotate the fractures
            if (AngleRotation == 0)
            {
                Frac_verts_host.insert(Frac_verts_host.end(), Frac_verts_host_i.begin(), Frac_verts_host_i.end());
                continue;
            }

            cuDFNsys::Quaternion<_DataType_> qua;
            qua = qua.DescribeRotation(rotationAxis, AngleRotation * M_PI / 180.0);

            for (uint j = 0; j < DSIZE; ++j)
            {
                for (uint k = 0; k < 4; ++k)
                {

                    cuDFNsys::Vector3<_DataType_> Vtex = Frac_verts_host_i[j].Verts3D[k];
                    Vtex.x -= Frac_verts_host_i[j].Center.x;
                    Vtex.y -= Frac_verts_host_i[j].Center.y;
                    Vtex.z -= Frac_verts_host_i[j].Center.z;

                    Vtex = qua.Rotate(Vtex);

                    Vtex.x += Frac_verts_host_i[j].Center.x;
                    Vtex.y += Frac_verts_host_i[j].Center.y;
                    Vtex.z += Frac_verts_host_i[j].Center.z;

                    Frac_verts_host_i[j].Verts3D[k] = Vtex;
                }

                cuDFNsys::Vector3<_DataType_> Vtex = Frac_verts_host_i[j].NormalVec;
                Vtex = qua.Rotate(Vtex);
                Frac_verts_host_i[j].NormalVec = Vtex;

                _DataType_ norm_f = sqrt(Frac_verts_host_i[j].NormalVec.x * Frac_verts_host_i[j].NormalVec.x +
                                         Frac_verts_host_i[j].NormalVec.y * Frac_verts_host_i[j].NormalVec.y +
                                         Frac_verts_host_i[j].NormalVec.z * Frac_verts_host_i[j].NormalVec.z);
                Frac_verts_host_i[j].NormalVec.x /= norm_f;
                Frac_verts_host_i[j].NormalVec.y /= norm_f;
                Frac_verts_host_i[j].NormalVec.z /= norm_f;

                if (Frac_verts_host_i[j].NormalVec.z < 0)
                {
                    Frac_verts_host_i[j].NormalVec.z *= -1;
                    Frac_verts_host_i[j].NormalVec.y *= -1;
                    Frac_verts_host_i[j].NormalVec.x *= -1;
                }

                Frac_verts_host_i[j].ConnectModelSurf[0] = 0,
                Frac_verts_host_i[j].ConnectModelSurf[1] = 0,
                Frac_verts_host_i[j].ConnectModelSurf[2] = 0,
                Frac_verts_host_i[j].ConnectModelSurf[3] = 0,
                Frac_verts_host_i[j].ConnectModelSurf[4] = 0,
                Frac_verts_host_i[j].ConnectModelSurf[5] = 0;
            }
            //------concatenate vectors
            Frac_verts_host.insert(Frac_verts_host.end(), Frac_verts_host_i.begin(), Frac_verts_host_i.end());
        }

        uint DSIZE = Frac_verts_host.size();

        thrust::device_vector<cuDFNsys::Fracture<_DataType_>> Frac_verts_device = Frac_verts_host;
        cuDFNsys::Fracture<_DataType_> *Frac_verts_device_ptr;
        Frac_verts_device_ptr = thrust::raw_pointer_cast(Frac_verts_device.data());

        cout << "dealing with all " << DSIZE << " fractures ..." << endl;
        cuDFNsys::FracturesChangeDomainSize<<<DSIZE / 256 + 1, 256>>>(Frac_verts_device_ptr,
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
        //for (size_t i = 0; i < ListClusters.size(); ++i)
        //    for (size_t j = 0; j < ListClusters[i].size(); ++j)
        //        cout << ListClusters[i][j] << (j == ListClusters[i].size() - 1 ? "\n" : ", ");

        cuDFNsys::IdentifyPercolationCluster<_DataType_> IdentiClu{ListClusters,
                                                                   Frac_verts_host, perco_dir,
                                                                   Percolation_cluster};
        cout << "DFN I finished" << endl;
        cuDFNsys::MatlabPlotDFN<_DataType_> As{"DFN_I.h5", "DFN_I.m",
                                               Frac_verts_host, Intersection_map, ListClusters,
                                               Percolation_cluster, false, true, true, true,
                                               L, perco_dir, true, "DFN_I"};

        cuDFNsys::OutputObjectData<_DataType_> lk;
        lk.OutputFractures("Fractures.h5", Frac_verts_host, L);

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
        //for (size_t i = 0; i < ListClusters.size(); ++i)
        //    for (size_t j = 0; j < ListClusters[i].size(); ++j)
        //        cout << ListClusters[i][j] << (j == ListClusters[i].size() - 1 ? "\n" : ", ");
        cout << "DFN II finished" << endl;
        cuDFNsys::MatlabPlotDFN<_DataType_> As2{"DFN_II.h5", "DFN_II.m",
                                                Frac_verts_host, Intersection_map, ListClusters,
                                                Percolation_cluster, true, true, true, true,
                                                L, perco_dir, true, "DFN_II"};
        //return 0;
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
            int NUMprior = Fracs_percol.size();

            bool ifRemoveDeadEnds = (atoi(argv[count_Argv_index + 1]) == 0 ? false : true);
            count_Argv_index++;

            cuDFNsys::RemoveDeadEndFrac<_DataType_> RDEF{Fracs_percol,
                                                         IntersectionPair_percol,
                                                         (size_t)perco_dir,
                                                         Frac_verts_host,
                                                         Intersection_map, ifRemoveDeadEnds};

            if (ifRemoveDeadEnds)
                cout << "remove " << NUMprior - Frac_verts_host.size() << " fractures\n";

            cout << "meshing ..." << endl;

            cuDFNsys::Mesh<_DataType_> mesh{Frac_verts_host, IntersectionPair_percol,
                                            &Fracs_percol,
                                            atof(argv[count_Argv_index + 1]),
                                            atof(argv[count_Argv_index + 2]),
                                            perco_dir,
                                            L};
            count_Argv_index += 2;

            lk.OutputMesh("mesh.h5", mesh, Fracs_percol);
            int i = 0;
            mesh.MatlabPlot("DFN_mesh_" + to_string(i + 1) + ".h5",
                            "DFN_mesh_" + to_string(i + 1) + ".m",
                            Frac_verts_host, L, true, true, true, "DFN_mesh_" + to_string(i + 1));

            cout << "MHFEM ing ..." << endl;

            cuDFNsys::MHFEM<_DataType_> fem{mesh, Frac_verts_host, 100, 20, perco_dir, L};
            lk.OutputMHFEM("mhfem.h5", fem);

            cout << "Fluxes: " << fem.QIn << ", ";
            cout << fem.QOut << ", Permeability: ";
            cout << fem.Permeability << endl;
            cout << "Error between the inlet and outlet fluxes: " << abs(fem.QIn - fem.QOut) / ((fem.QOut + fem.QIn) * 0.5) * 100.0 << "%\n";
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

            cout << "Particle transport ing ...\n";
            return 0;
            // cuDFNsys::ParticleTransport<_DataType_> p{atoi(argv[13]),             // number of particle
            //                                           atoi(argv[14]),             // number of time steps
            //                                           (_DataType_)atof(argv[15]), // delta T
            //                                           (_DataType_)atof(argv[16]), // molecular diffusion
            //                                           Frac_verts_host, mesh, fem, (uint)perco_dir, -0.5f * L,
            //                                           "Particle_tracking", "Flux-weighted"};
            // p.MatlabPlot("MHFEM_" + to_string(i + 1) + ".h5", "particle.m", mesh, fem, L);
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