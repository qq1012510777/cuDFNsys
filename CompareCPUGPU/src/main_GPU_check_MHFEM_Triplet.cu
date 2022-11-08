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

int main(int argc, char *argv[])
{

    try
    {

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

            //int Nproc = atoi(argv[1]);
            uint MonteCarloTimes = atoi(argv[1]);

            uint2 Dims = make_uint2(1, MonteCarloTimes);

            string filename = "GPU_check_MHFEM_Triplet";

            cuDFNsys::HDF5API hd5;

            string hdfilename = filename + "/" + "MHFEM_Triplet_countTime.h5";

            for (size_t k = 1; k <= 20; k++)
            {
                if (k == 1)
                {
                    int us = remove(filename.c_str());

                    us = mkdir(filename.c_str(), 0777);

                    if (us == -1)
                        throw cuDFNsys::ExceptionsPause("Cannot create file " + filename);

                    hd5.NewFile(hdfilename);

                    hd5.AddDataset(hdfilename, "N",
                                   "MonteCarloTimes", &MonteCarloTimes, make_uint2(1, 1));
                }

                _DataType_ minGrid = 21 - k * 1.0;
                _DataType_ maxGrid = minGrid + 1;

                hd5.AddDataset(hdfilename, "Step_" + cuDFNsys::ToStringWithWidth(k, 5),
                               "minGrid", &minGrid, make_uint2(1, 1));
                hd5.AddDataset(hdfilename, "Step_" + cuDFNsys::ToStringWithWidth(k, 5),
                               "maxGrid", &maxGrid, make_uint2(1, 1));

                vector<double> TripletTime(MonteCarloTimes);
                vector<double> Q_in(MonteCarloTimes);
                vector<double> Q_out(MonteCarloTimes);
                vector<double> NUMeles(MonteCarloTimes);

                for (uint op = 0; op < MonteCarloTimes; ++op)
                {
                    cout << "meshing ..." << k << "/" << op + 1 << endl;

                    cuDFNsys::Mesh<_DataType_> mesh{Frac_verts_host, IntersectionPair_percol,
                                                    &Fracs_percol, minGrid, maxGrid, perco_dir, L};

                    cout << "MHFEM ing ..." << k << "/" << op + 1 << endl;

                    cuDFNsys::MHFEM<_DataType_> fem{mesh, Frac_verts_host, 100, 20, perco_dir, L};
                    cout << "Fluxes: " << fem.QIn << ", ";
                    cout << fem.QOut << ", Permeability: ";
                    cout << fem.Permeability << endl;
                    if (fem.QError > 1 || isnan(fem.Permeability) == 1)
                        cout << "\e[1;32mFound large error or isnan, the error: " << fem.QError << ", the permeability: " << fem.Permeability << "\e[0m\n";

                    TripletTime[op] = fem.TripletTime;
                    Q_in[op] = fem.QIn;
                    Q_out[op] = fem.QOut;
                    NUMeles[op] = mesh.Element3D.size();
                }

                hd5.AddDataset(hdfilename, "Step_" + cuDFNsys::ToStringWithWidth(k, 5),
                               "TripletTime", TripletTime.data(), Dims);
                hd5.AddDataset(hdfilename, "Step_" + cuDFNsys::ToStringWithWidth(k, 5),
                               "Q_in", Q_in.data(), Dims);
                hd5.AddDataset(hdfilename, "Step_" + cuDFNsys::ToStringWithWidth(k, 5),
                               "Q_out", Q_out.data(), Dims);
                hd5.AddDataset(hdfilename, "Step_" + cuDFNsys::ToStringWithWidth(k, 5),
                               "NUMeles", NUMeles.data(), Dims);
            }
        }
        //cudaDeviceReset();
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