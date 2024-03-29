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
        // to validate dispersion in two intersected fractures / or two non-intersected fractures
        // Re = 0.01 LengthScale = 30
        // NUMstep = 50
        // Factor_mean_time_in_grid = 1e5
        int NumParticles = atoi(argv[1]);
        double Factor_mean_time_in_grid = atof(argv[2]);
        double LengthScale_Over_Pe = atof(argv[3]) / atof(argv[4]);
        uint NUMstep = atoi(argv[5]);
        bool IfSetOtherInjectPlane = (atoi(argv[6]) == 0 ? false : true);
        _DataType_ InjectPlane = atof(argv[7]);

        double istart = cuDFNsys::CPUSecond();

        int dev = 0;
        GPUErrCheck(cudaSetDevice(dev));

        int DSIZE = 0;
        _DataType_ L = 0;
        _DataType_ minGrid = 0;
        _DataType_ maxGrid = 0;

        L = 30;
        DSIZE = 800; //------------------------------------------
        minGrid = 1;
        maxGrid = 3;

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
        if (DSIZE == 1)
        {
            cuDFNsys::FracturesCrossedVertical<_DataType_><<<DSIZE / 256 + 1, 256>>>(Frac_verts_device_ptr,
                                                                                     (unsigned long)t,
                                                                                     DSIZE,
                                                                                     L);
            // cuDFNsys::FracturesIncomplete<_DataType_><<<DSIZE / 256 + 1, 256>>>(Frac_verts_device_ptr,
            //                                                                     (unsigned long)t,
            //                                                                     DSIZE,
            //                                                                     L);
            // cuDFNsys::FractureTwoIntersectOrNot<_DataType_><<<DSIZE / 256 + 1, 256>>>(Frac_verts_device_ptr,
            //                                                                           (unsigned long)t,
            //                                                                           DSIZE,
            //                                                                           L,
            //                                                                           true);
        }
        else if (DSIZE == 2)
            cuDFNsys::FractureTwoIntersectOrNot<_DataType_><<<DSIZE / 256 + 1, 256>>>(Frac_verts_device_ptr,
                                                                                      (unsigned long)t,
                                                                                      DSIZE,
                                                                                      L,
                                                                                      true);
        else
            cuDFNsys::FracturesParallel<_DataType_><<<DSIZE / 256 + 1, 256>>>(Frac_verts_device_ptr,
                                                                              DSIZE,
                                                                              (unsigned long)(t * t),
                                                                              L);

        cudaDeviceSynchronize();

        Frac_verts_host = Frac_verts_device;

        thrust::host_vector<cuDFNsys::Fracture<_DataType_>> Frac_verts_host_II;

        for (int i = 0; i < Frac_verts_host.size(); ++i)
            if (Frac_verts_host[i].ConnectModelSurf[4] && Frac_verts_host[i].ConnectModelSurf[5])
                Frac_verts_host_II.push_back(Frac_verts_host[i]);

        Frac_verts_host = Frac_verts_host_II;

        DSIZE = Frac_verts_host.size();
        cout << "DSIZE: " << DSIZE << endl;
        std::map<pair<size_t, size_t>, pair<cuDFNsys::Vector3<_DataType_>, cuDFNsys::Vector3<_DataType_>>> Intersection_map;
        // cout << "Frac_verts_host.size(): " << Frac_verts_host.size() << endl;
        // cout << Frac_verts_host[0].Verts3D[0].x << ", " << Frac_verts_host[0].Verts3D[0].y << ", " << Frac_verts_host[0].Verts3D[0].z << endl;
        // cout << Frac_verts_host[0].Verts3D[1].x << ", " << Frac_verts_host[0].Verts3D[1].y << ", " << Frac_verts_host[1].Verts3D[0].z << endl;
        // cout << Frac_verts_host[0].Verts3D[2].x << ", " << Frac_verts_host[0].Verts3D[2].y << ", " << Frac_verts_host[2].Verts3D[0].z << endl;
        // cout << Frac_verts_host[0].Verts3D[3].x << ", " << Frac_verts_host[0].Verts3D[3].y << ", " << Frac_verts_host[3].Verts3D[0].z << endl;

        if (DSIZE > 1)
        {
            cout << "identifying intersections with truncated fractures" << endl;
            cuDFNsys::IdentifyIntersection<_DataType_> identifyInters{Frac_verts_host.size(),
                                                                      Frac_verts_device_ptr,
                                                                      true,
                                                                      Intersection_map};
        }

        std::vector<std::vector<size_t>> ListClusters;
        std::vector<size_t> Percolation_cluster;

        if (DSIZE == 1)
        {
            ListClusters.resize(1);
            Percolation_cluster.resize(1);
            ListClusters[0] = std::vector<size_t>(0);
            Percolation_cluster[0] = 0;
        }
        else
        {
            cout << "identifying cluster with truncated fractures" << endl;
            cuDFNsys::Graph<_DataType_> G2{(size_t)DSIZE, Intersection_map};
            G2.UseDFS(ListClusters);
            cuDFNsys::IdentifyPercolationCluster<_DataType_> IdentiClu2{ListClusters,
                                                                        Frac_verts_host, perco_dir,
                                                                        Percolation_cluster};
            // for (auto e : Percolation_cluster)
            // cout << e << endl;
        }

        cuDFNsys::MatlabPlotDFN<_DataType_> As2{"DFN_II.h5",
                                                "DFN_II.m",
                                                Frac_verts_host,
                                                Intersection_map,
                                                ListClusters,
                                                Percolation_cluster,
                                                true, true, true, true,
                                                L, perco_dir, true, "DFN_II"};

        Frac_verts_device.clear();
        Frac_verts_device.shrink_to_fit();
        //-----------
        //exit(0);

        cuDFNsys::OutputObjectData<_DataType_> lk;
        if (Percolation_cluster.size() > 0)
        {

            double istart_1 = cuDFNsys::CPUSecond();

            std::vector<size_t> Fracs_percol;
            std::vector<pair<int, int>> IntersectionPair_percol;

            if (DSIZE == 2)
            {
                cuDFNsys::GetAllPercolatingFractures GetPer{Percolation_cluster,
                                                            ListClusters,
                                                            Fracs_percol};

                cuDFNsys::RemoveDeadEndFrac<_DataType_> RDEF{Fracs_percol,
                                                             IntersectionPair_percol,
                                                             (size_t)perco_dir,
                                                             Frac_verts_host,
                                                             Intersection_map,
                                                             false};
                // cout << Fracs_percol.size() << endl;
                // cout << IntersectionPair_percol.size() << endl;
                IntersectionPair_percol.resize(1);
                IntersectionPair_percol[0] = std::make_pair(0, 1);
            }
            else if (DSIZE == 1)
            {
                Fracs_percol = Percolation_cluster;
                IntersectionPair_percol.resize(0);
            }
            else
            {
                cuDFNsys::GetAllPercolatingFractures GetPer{Percolation_cluster,
                                                            ListClusters,
                                                            Fracs_percol};

                cuDFNsys::RemoveDeadEndFrac<_DataType_> RDEF{Fracs_percol,
                                                             IntersectionPair_percol,
                                                             (size_t)perco_dir,
                                                             Frac_verts_host,
                                                             Intersection_map,
                                                             false};
                //Fracs_percol.resize(Frac_verts_host.size());

                IntersectionPair_percol.resize(0);
                // for (int e = 0; e < IntersectionPair_percol.size(); ++e)
                // {
                //     //Fracs_percol[e] = e;
                //     IntersectionPair_percol[e] = std::make_pair(Fracs_percol[e], Fracs_percol[e]);
                // }
            }

            cout << "meshing ..." << endl;

            cuDFNsys::Mesh<_DataType_> mesh{Frac_verts_host, IntersectionPair_percol,
                                            &Fracs_percol, minGrid, maxGrid, perco_dir, L};

            lk.OutputMesh("mesh.h5", mesh, Fracs_percol);

            //lk.OutputFractures("FracturesII.h5", Frac_verts_host, L, make_double3(1, 1, 1));

            cout << "MHFEM ing ..." << endl;

            cuDFNsys::MHFEM<_DataType_> fem{mesh, Frac_verts_host, L, 0.0, perco_dir, L};
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
            double mean_grid_area = mesh.MatlabPlot("DFN_mesh_" + to_string(i + 1) + ".h5",
                                                    "DFN_mesh.m",
                                                    Frac_verts_host, L, true, true, true, "DFN_mesh_" + to_string(i + 1));
            double2 TGH = fem.MatlabPlot("MHFEM_" + to_string(i + 1) + ".h5",
                                         "MHFEM_" + to_string(i + 1) + ".m",
                                         Frac_verts_host, mesh, L, true, "MHFEM_" + to_string(i + 1));
            //---------------
            cout << "The mean area of all elements is " << mean_grid_area << endl;
            double meanV = TGH.x;
            double maxV = TGH.y;
            cout << "The maximum velocity of all elements is " << maxV << endl;
            cout << "The mean velocity of all elements is " << meanV << endl;

            double DiffusionLocal = LengthScale_Over_Pe * meanV;
            double meanTime = pow(mean_grid_area, 0.5) / maxV;
            double DeltaT = meanTime / Factor_mean_time_in_grid;

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

            cuDFNsys::ParticleTransport<_DataType_> p{(int)NUMstep, // number of time step
                                                      Frac_verts_host,
                                                      mesh,
                                                      fem,
                                                      (uint)perco_dir,
                                                      -0.5f * L,
                                                      NumParticles, // num of particle
                                                      DeltaT,       // delta_T_ii
                                                      DiffusionLocal,
                                                      "Particle_tracking",
                                                      "Resident",
                                                      "OutputAll",
                                                      false,
                                                      1,
                                                      false,
                                                      1000,
                                                      true,
                                                      IfSetOtherInjectPlane,
                                                      InjectPlane, false};

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