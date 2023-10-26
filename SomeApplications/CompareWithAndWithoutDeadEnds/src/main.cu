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
// NAME:        CompareWithAndWithoutDeadEnds
// DESCRIPTION: Compare permeability With And Without DeadEnds
// AUTHOR:      Tingchang YIN
// DATE:        23/11/2022
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
    int iDev = 0;
    GPUErrCheck(cudaSetDevice(iDev));

    double iStart = cuDFNsys::CPUSecond();

    uint IfFlowModel = 0;
    IfFlowModel = (uint)(atoi(argv[1]));

    uint ModeSizeDistr = (uint)(atoi(argv[2]));
    cuDFNsys::Vector4<_DataType_> SizeParaDis = cuDFNsys::MakeVector4(atof(argv[3]),
                                                                      atof(argv[4]),
                                                                      atof(argv[5]),
                                                                      atof(argv[6]));

    _DataType_ kappa = atof(argv[7]);
    _DataType_ conductivity_powerlaw_e = atof(argv[8]);

    _DataType_ L = atof(argv[9]);
    int perco_dir = atoi(argv[10]);

    int Init_NUM_Frac = atoi(argv[11]);

    int Frac_increment = atoi(argv[12]);

    int LOOP_times = atoi(argv[13]);
    int inti_LOOP_times = 0;

    _DataType_ minGrid = atof(argv[14]);
    _DataType_ maxGrid = atof(argv[15]);
    _DataType_ InletP = atof(argv[16]);
    _DataType_ OutletP = atof(argv[17]);

    int MODEL_NO = atoi(argv[18]);
    int MC_NO = atoi(argv[19]);

    _DataType_ Gamma_constant = 5.0e-4;

    if (argv[20] != NULL)
    {
        Gamma_constant = atof(argv[20]);
        cout << "\n\n********** cuDFNsys receives non-default gamma=" << Gamma_constant << " constant; b = gamma * R ^ beta**********\n\n";
    }
    else
        cout << "\n\n********** cuDFNsys uses default gamma constant; b = gamma * R ^ beta**********\n\n";

    _DataType_ P33_total_A[1] = {0};
    _DataType_ P33_connected_A[1] = {0};
    _DataType_ Ratio_of_P33_A[1] = {0};
    _DataType_ P33_largest_cluster_A[1] = {0};

    _DataType_ P32_total_A[1] = {0};
    _DataType_ P32_connected_A[1] = {0};
    _DataType_ Ratio_of_P32_A[1] = {0};
    _DataType_ P32_largest_cluster_A[1] = {0};

    _DataType_ P30_A[1] = {0};
    _DataType_ P30_connected_A[1] = {0};
    _DataType_ Ratio_of_P30_A[1] = {0};
    _DataType_ P30_largest_cluster_A[1] = {0};

    _DataType_ Percolation_probability_A[1] = {0};
    _DataType_ n_I_A[1] = {0};

    _DataType_ Permeability_A[1] = {0};
    _DataType_ Q_error_A[1] = {0};

    _DataType_ NumFractureRemoved[1] = {-1};

    cuDFNsys::HDF5API h5out;
    string filename = "Datafile_Model_" + cuDFNsys::ToStringWithWidth(MODEL_NO, 3) +
                      "_MC_" + cuDFNsys::ToStringWithWidth(MC_NO, 5) + ".h5";

    uint2 dim_e = make_uint2(1, 1);

    if (!h5out.IfH5FileExist(filename))
    {
        h5out.NewFile(filename);
        _DataType_ domainsize[1] = {L};
        h5out.AddDataset(filename, "N", "Domain_size", domainsize, dim_e);
    }
    else
    {
        try
        {
            vector<double> Looptimes_k =
                h5out.ReadDataset<double>(filename, "N", "Loop_times");

            if (Looptimes_k[0] >= LOOP_times)
            {
                cout << "This simulation has been finished!\n";
                return 0;
            }
            else
            {
                inti_LOOP_times = (int)Looptimes_k[0];
            }
        }
        catch (...)
        {
            // no dataset "Loop_times" existing means that I did not start loop
            //
            inti_LOOP_times = 0;
        }
    }

    for (int i = inti_LOOP_times; i < LOOP_times; ++i)
    {
    ReGen_111:;
        try
        {
            cout << "LOOP " << i << endl;
            double iStart_DFN = cuDFNsys::CPUSecond();

            int DSIZE = Init_NUM_Frac + (i + 1) * Frac_increment;

            thrust::host_vector<cuDFNsys::Fracture<_DataType_>> Frac_verts_host(DSIZE);
            thrust::device_vector<cuDFNsys::Fracture<_DataType_>> Frac_verts_device(DSIZE);
            cuDFNsys::Fracture<_DataType_> *Frac_verts_device_ptr = thrust::raw_pointer_cast(Frac_verts_device.data());
            //cout << "\tAllocated memory to vectors of cuDFNsys::Fracture<_DataType_>\n";
            cuDFNsys::Warmup<<<DSIZE / 256 + 1, 256 /*  1, 2*/>>>();
            cudaDeviceSynchronize();
            time_t t;
            time(&t);
            //cout << "conductivity_powerlaw_e: " << conductivity_powerlaw_e << endl;
            cuDFNsys::Fractures<_DataType_><<<DSIZE / 256 + 1, 256 /*  1, 2*/>>>(Frac_verts_device_ptr,
                                                                                 (unsigned long)t,
                                                                                 DSIZE,
                                                                                 L,
                                                                                 ModeSizeDistr,
                                                                                 SizeParaDis,
                                                                                 kappa,
                                                                                 conductivity_powerlaw_e, Gamma_constant);
            //-----------
            cudaDeviceSynchronize();
            //cout << "\tGenerated fractures\n";
            Frac_verts_host = Frac_verts_device;

            std::map<pair<size_t, size_t>, pair<cuDFNsys::Vector3<_DataType_>, cuDFNsys::Vector3<_DataType_>>> Intersection_map;
            try
            {
                cuDFNsys::IdentifyIntersection<_DataType_> identifyInters{Frac_verts_host.size(),
                                                                          Frac_verts_device_ptr,
                                                                          false,
                                                                          Intersection_map};
            }
            catch (thrust::system::system_error &e)
            {
                cout << "\033[1;32m\tThrew an exception: " << e.what() << "\n\tLet's use CPU finish identification of intersections\033[0m" << endl;
                cudaDeviceReset();
                Intersection_map.clear();
                cuDFNsys::IdentifyIntersection<_DataType_> identifyInters3{Frac_verts_host,
                                                                           false,
                                                                           Intersection_map};
            };
            //cout << "\tIdentified intersection of DFN with complete factures\n";
            std::vector<std::vector<size_t>> ListClusters;
            std::vector<size_t> Percolation_cluster;
            cuDFNsys::Graph<_DataType_> G{(size_t)DSIZE, Intersection_map};
            G.UseDFS(ListClusters);
            cuDFNsys::IdentifyPercolationCluster<_DataType_> IdentiClu{ListClusters,
                                                                       Frac_verts_host, perco_dir,
                                                                       Percolation_cluster};

            cuDFNsys::GetStatistics<_DataType_>(Frac_verts_host,
                                                Intersection_map,
                                                ListClusters,
                                                Percolation_cluster,
                                                L,
                                                P33_total_A[0],
                                                P33_connected_A[0],
                                                Ratio_of_P33_A[0],
                                                P33_largest_cluster_A[0],
                                                P32_total_A[0],
                                                P32_connected_A[0],
                                                Ratio_of_P32_A[0],
                                                P32_largest_cluster_A[0],
                                                P30_A[0],
                                                P30_connected_A[0],
                                                Ratio_of_P30_A[0],
                                                P30_largest_cluster_A[0],
                                                Percolation_probability_A[0],
                                                n_I_A[0]);

            Intersection_map.clear();
            ListClusters.clear();
            Percolation_cluster.clear();
            try
            {
                cuDFNsys::IdentifyIntersection<_DataType_> identifyInters2{Frac_verts_host.size(),
                                                                           Frac_verts_device_ptr, true,
                                                                           Intersection_map};
            }
            catch (thrust::system::system_error &e)
            {
                cout << "\033[1;32m\tThrew an exception: " << e.what() << "\n\tLet's use CPU finish identification of intersections\033[0m" << endl;
                cudaDeviceReset();
                Intersection_map.clear();
                cuDFNsys::IdentifyIntersection<_DataType_> identifyInters3{Frac_verts_host,
                                                                           true,
                                                                           Intersection_map};
            };
            //cout << "\tIdentified intersection of DFN with truncated factures\n";
            cuDFNsys::Graph<_DataType_> G2{(size_t)DSIZE, Intersection_map};
            G2.UseDFS(ListClusters);
            cuDFNsys::IdentifyPercolationCluster<_DataType_> IdentiClu2{ListClusters,
                                                                        Frac_verts_host, perco_dir,
                                                                        Percolation_cluster};

            double iElaps_DFN = cuDFNsys::CPUSecond() - iStart_DFN;
            cout << "\tDFN generated. Running time: " << iElaps_DFN << " sec\n";

            Frac_verts_device.clear();
            Frac_verts_device.shrink_to_fit();
            //cudaDeviceReset();

            //Connections.clear();
            //Connections.shrink_to_fit();
            thrust::host_vector<cuDFNsys::Fracture<_DataType_>> Frac_verts_host_back_up = Frac_verts_host;

            //-------------------------
            bool IfRemoveDeadEnds = false;
            if (IfFlowModel == 1)
            {

                if (Percolation_cluster.size() > 0)
                {
                    double iStart_mesh = cuDFNsys::CPUSecond();
                    std::vector<size_t> Fracs_percol;
                    cuDFNsys::GetAllPercolatingFractures GetPer{Percolation_cluster,
                                                                ListClusters,
                                                                Fracs_percol};

                    std::vector<pair<int, int>> IntersectionPair_percol;

                    int NUMFrac_perco_prior_remove = Fracs_percol.size();

                    cuDFNsys::RemoveDeadEndFrac<_DataType_> RDEF{Fracs_percol,
                                                                 IntersectionPair_percol,
                                                                 (size_t)perco_dir,
                                                                 Frac_verts_host,
                                                                 Intersection_map,
                                                                 IfRemoveDeadEnds};
                    if (IfRemoveDeadEnds)
                        NumFractureRemoved[0] = NUMFrac_perco_prior_remove - Frac_verts_host.size();

                    if (!IfRemoveDeadEnds) // do not remove
                    {
                        if (NUMFrac_perco_prior_remove != Frac_verts_host.size())
                        {
                            throw cuDFNsys::ExceptionsPause("The fractures in percolating clusters should not be removed!");
                        }
                    }

                    cout << "\tMeshing ... with Dead-ends ..." << endl;
                    //cudaDeviceReset();
                    cuDFNsys::Mesh<_DataType_> mesh{Frac_verts_host, IntersectionPair_percol,
                                                    &Fracs_percol, minGrid, maxGrid, perco_dir, L};

                    double iElaps_mesh = cuDFNsys::CPUSecond() - iStart_mesh;
                    cout << "\tMesh finished. Running time: " << iElaps_mesh << " sec\n";

                    cout << "\tMHFEM ing ..." << endl;
                    double iStart_mhfem = cuDFNsys::CPUSecond();
                    cuDFNsys::MHFEM<_DataType_> fem{mesh, Frac_verts_host, InletP, OutletP, perco_dir, L};

                    //cout << "\tFluxes: " << fem.Q_in << ", " << fem.Q_out << ", Permeability: " << fem.Permeability << endl;
                    if (fem.QError > 1 || isnan(fem.Permeability) == 1)
                    {
                        string AS = "Found large error or isnan, the error: " + to_string(fem.QError) + ", the permeability: " + to_string(fem.Permeability);
                        throw cuDFNsys::ExceptionsIgnore(AS);
                    }

                    Permeability_A[0] = fem.Permeability;
                    Q_error_A[0] = fem.QError;

                    double iElaps_mhfem = cuDFNsys::CPUSecond() - iStart_mhfem;
                    cout << "\tMHFEM finished. Running time: " << iElaps_mhfem << " sec\n";
                    //---------------------
                }
            }
            cout << "\tOutputing data... with dead-ends ...\n";
            string groupname = "group_" + cuDFNsys::ToStringWithWidth(i + 1, 3) + "_withDeadEnds";
            vector<string> datasetname = {
                "P33_total",
                "P33_connected_z",
                "Ratio_of_P33_z",
                "P33_largest_cluster",
                "P32_total",
                "P32_connected_z",
                "Ratio_of_P32_z",
                "P32_largest_cluster",
                "P30",
                "P30_connected_z",
                "Ratio_of_P30_z",
                "P30_largest_cluster",
                "Percolation_probability_z",
                "n_I",
                "Permeability_z",
                "Q_error_z"};
            vector<_DataType_ *> data_input = {P33_total_A,
                                               P33_connected_A,
                                               Ratio_of_P33_A,
                                               P33_largest_cluster_A,
                                               P32_total_A,
                                               P32_connected_A,
                                               Ratio_of_P32_A,
                                               P32_largest_cluster_A,
                                               P30_A,
                                               P30_connected_A,
                                               Ratio_of_P30_A,
                                               P30_largest_cluster_A,
                                               Percolation_probability_A,
                                               n_I_A,
                                               Permeability_A,
                                               Q_error_A};
            uint2 dim_ki;
            dim_ki.x = 1;
            dim_ki.y = 1;
            vector<uint2> dim_ss(data_input.size(), dim_ki);
            h5out.AddDatasetsWithOneGroup(filename, groupname,
                                          datasetname, data_input, dim_ss);

            if (IfRemoveDeadEnds)
                h5out.AddDataset(filename, groupname, "NumFractureRemoved", NumFractureRemoved, dim_e);

            data_input.clear();
            data_input.shrink_to_fit();

            //-------------------------- remove dead ends
            //-------------------------- remove dead ends
            //-------------------------- remove dead ends
            Frac_verts_host.clear();
            Frac_verts_host.shrink_to_fit();

            Frac_verts_host = Frac_verts_host_back_up;

            Frac_verts_host_back_up.clear();
            Frac_verts_host_back_up.shrink_to_fit();

            IfRemoveDeadEnds = true;
            if (IfFlowModel == 1)
            {

                if (Percolation_cluster.size() > 0)
                {
                    double iStart_mesh = cuDFNsys::CPUSecond();
                    std::vector<size_t> Fracs_percol;
                    cuDFNsys::GetAllPercolatingFractures GetPer{Percolation_cluster,
                                                                ListClusters,
                                                                Fracs_percol};

                    std::vector<pair<int, int>> IntersectionPair_percol;

                    int NUMFrac_perco_prior_remove = Fracs_percol.size();

                    cuDFNsys::RemoveDeadEndFrac<_DataType_> RDEF{Fracs_percol,
                                                                 IntersectionPair_percol,
                                                                 (size_t)perco_dir,
                                                                 Frac_verts_host,
                                                                 Intersection_map,
                                                                 IfRemoveDeadEnds};
                    if (IfRemoveDeadEnds)
                        NumFractureRemoved[0] = NUMFrac_perco_prior_remove - Frac_verts_host.size();

                    if (!IfRemoveDeadEnds) // do not remove
                    {
                        if (NUMFrac_perco_prior_remove != Frac_verts_host.size())
                        {
                            throw cuDFNsys::ExceptionsPause("The fractures in percolating clusters should not be removed!");
                        }
                    }

                    cout << "\tMeshing ... rm Dead-ends ..." << endl;
                    //cudaDeviceReset();
                    cuDFNsys::Mesh<_DataType_> mesh{Frac_verts_host, IntersectionPair_percol,
                                                    &Fracs_percol, minGrid, maxGrid, perco_dir, L};

                    double iElaps_mesh = cuDFNsys::CPUSecond() - iStart_mesh;
                    cout << "\tMesh finished. Running time: " << iElaps_mesh << " sec\n";

                    cout << "\tMHFEM ing ..." << endl;
                    double iStart_mhfem = cuDFNsys::CPUSecond();
                    cuDFNsys::MHFEM<_DataType_> fem{mesh, Frac_verts_host, InletP, OutletP, perco_dir, L};

                    //cout << "\tFluxes: " << fem.Q_in << ", " << fem.Q_out << ", Permeability: " << fem.Permeability << endl;
                    if (fem.QError > 1 || isnan(fem.Permeability) == 1)
                    {
                        string AS = "Found large error or isnan, the error: " + to_string(fem.QError) + ", the permeability: " + to_string(fem.Permeability);
                        throw cuDFNsys::ExceptionsIgnore(AS);
                    }

                    Permeability_A[0] = fem.Permeability;
                    Q_error_A[0] = fem.QError;

                    double iElaps_mhfem = cuDFNsys::CPUSecond() - iStart_mhfem;
                    cout << "\tMHFEM finished. Running time: " << iElaps_mhfem << " sec\n";
                    //---------------------
                }
            }
            cout << "\tOutputing data... rm dead-ends ...\n";
            groupname = "group_" + cuDFNsys::ToStringWithWidth(i + 1, 3) + "_rmDeadEnds";
            vector<_DataType_ *> data_inputII = {P33_total_A,
                                                 P33_connected_A,
                                                 Ratio_of_P33_A,
                                                 P33_largest_cluster_A,
                                                 P32_total_A,
                                                 P32_connected_A,
                                                 Ratio_of_P32_A,
                                                 P32_largest_cluster_A,
                                                 P30_A,
                                                 P30_connected_A,
                                                 Ratio_of_P30_A,
                                                 P30_largest_cluster_A,
                                                 Percolation_probability_A,
                                                 n_I_A,
                                                 Permeability_A,
                                                 Q_error_A};
            h5out.AddDatasetsWithOneGroup(filename, groupname,
                                          datasetname, data_inputII, dim_ss);

            //----------------------------
            //----------------------------
            //----------------------------

            _DataType_ i_p[1] = {(_DataType_)i + 1};
            if (i == 0)
                h5out.AddDataset(filename, "N", "Loop_times", i_p, dim_e);
            else
                h5out.OverWrite(filename, "N", "Loop_times", i_p, dim_e);
        }
        catch (cuDFNsys::ExceptionsIgnore &e)
        {
            cout << e.what() << endl;
            goto ReGen_111;
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
    };

    cout << "Loop finished!\n";
    double iElaps = cuDFNsys::CPUSecond() - iStart;
    cout << "Running times of this loop is: " << iElaps << " sec\n";
    return 0;
};