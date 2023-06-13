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
// NAME:        PercolationAnisotropicCases
// DESCRIPTION: PercolationAnisotropicCases
// AUTHOR:      Tingchang YIN
// DATE:        26/05/2023
// ====================================================

#include "cuDFNsys.cuh"
#include <fstream>
#include <iostream>
#include <limits.h>
#include <unistd.h>

#ifdef USE_DOUBLES
typedef double _DataType_;
#else
typedef float _DataType_;
#endif
// _DataType_ here is actually double

int main(int argc, char *argv[])
{

    try
    {
        int iDev = 0;
        GPUErrCheck(cudaSetDevice(iDev));

        bool IfFlowSimulation = atoi(argv[1]) == 0 ? false : true;

        uint ModeOfFractureSizeDistribution = atoi(argv[2]);
        if (ModeOfFractureSizeDistribution > 3)
            throw cuDFNsys::ExceptionsPause("Wrong mode of fracture size distributions\n");

        cuDFNsys::Vector4<_DataType_> SizeDistributionParameters = cuDFNsys::MakeVector4(atof(argv[3]),
                                                                                         atof(argv[4]),
                                                                                         atof(argv[5]),
                                                                                         atof(argv[6]));
        _DataType_ kappa = atof(argv[7]);

        _DataType_ conductivity_beta = atof(argv[8]);

        _DataType_ conductivity_gamma = atof(argv[9]);

        if (conductivity_gamma == 5e-4)
            throw cuDFNsys::ExceptionsPause("Don't set conductivity_gamma = 5e-4\n");

        _DataType_ L = atof(argv[10]);

        double3 DomainDimensionRatio = make_double3(atof(argv[11]),
                                                    atof(argv[12]),
                                                    atof(argv[13]));

        int perco_dir = atoi(argv[14]);

        int InitNumFractures = atoi(argv[15]);

        int IncrementFractures = atoi(argv[16]);

        int LOOP_times = atoi(argv[17]);
        int inti_LOOP_times = 0;

        _DataType_ minGrid = atof(argv[18]);
        _DataType_ maxGrid = atof(argv[19]);

        int MODEL_NO = atoi(argv[20]);
        int MC_NO = atoi(argv[21]);

        cout << "//-----------------------input paramters----------------------------\n";
        cout << "IfFlowSimulation: " << (IfFlowSimulation ? "true" : "false") << endl;
        cout << "ModeOfFractureSizeDistribution: " << ModeOfFractureSizeDistribution << endl;
        cout << "SizeDistributionParameters: " << SizeDistributionParameters.x << ", "
             << SizeDistributionParameters.y << ", "
             << SizeDistributionParameters.z << ", "
             << SizeDistributionParameters.w << endl;
        cout << "kappa: " << kappa << endl;
        cout << "conductivity_beta: " << conductivity_beta << endl;
        cout << "conductivity_gamma: " << conductivity_gamma << endl;
        cout << "L: " << L << endl;
        cout << "DomainDimensionRatio: " << DomainDimensionRatio.x << ", "
             << DomainDimensionRatio.y << ", "
             << DomainDimensionRatio.z << endl;
        cout << "perco_dir: " << perco_dir << endl;

        cout << "InitNumFractures: " << InitNumFractures << endl;
        cout << "IncrementFractures: " << IncrementFractures << endl;
        cout << "LOOP_times: " << LOOP_times << endl;
        cout << "minGrid: " << minGrid << endl;
        cout << "maxGrid: " << maxGrid << endl;
        cout << "MODEL_NO: " << MODEL_NO << endl;
        cout << "MC_NO: " << MC_NO << endl;
        cout << "//-------------------------------------------------------------------\n\n";
        //----------------------------------
        //----------------------------------
        //----------------------------------

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

        cuDFNsys::HDF5API h5out;
        string filename = "Datafile_Model_" + cuDFNsys::ToStringWithWidth(MODEL_NO, 3) +
                          "_MC_" + cuDFNsys::ToStringWithWidth(MC_NO, 5) + ".h5";

        uint2 dim_e = make_uint2(1, 1);

        if (!h5out.IfH5FileExist(filename))
        {
            h5out.NewFile(filename);
            _DataType_ AspectRatio[3] = {DomainDimensionRatio.x,
                                         DomainDimensionRatio.y,
                                         DomainDimensionRatio.z};
            h5out.AddDataset(filename, "N", "Lm", &L, dim_e);
            h5out.AddDataset(filename, "N", "PercolationDirection",
                             &perco_dir, dim_e);
            h5out.AddDataset(filename, "N", "DomainDimensionRatio",
                             AspectRatio, make_uint2(3, 1));
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

        //------------------------
        //------------------------
        //------------------------

        for (int i = inti_LOOP_times; i < LOOP_times; ++i)
        {
        ReGen_111:;

            try
            {
                cout << "LOOP " << i << endl;

                double iStart_DFN = cuDFNsys::CPUSecond();

                int DSIZE = InitNumFractures + (i + 1) * IncrementFractures;

                thrust::host_vector<cuDFNsys::Fracture<_DataType_>> Frac_verts_host(DSIZE);
                thrust::device_vector<cuDFNsys::Fracture<_DataType_>> Frac_verts_device(DSIZE);
                cuDFNsys::Fracture<_DataType_> *Frac_verts_device_ptr = thrust::raw_pointer_cast(Frac_verts_device.data());
                //cout << "\tAllocated memory to vectors of cuDFNsys::Fracture<_DataType_>\n";
                cuDFNsys::Warmup<<<DSIZE / 256 + 1, 256 /*  1, 2*/>>>();
                cudaDeviceSynchronize();
                time_t t;
                time(&t);

                Eigen::MatrixXd Ter = Eigen::MatrixXd::Random(1, 1);
                unsigned long TYR = (unsigned long)t + (unsigned long)ceil(abs(Ter(0, 0)) * ((unsigned long)t * 1.0));

                cuDFNsys::Fractures<_DataType_><<<DSIZE / 256 + 1, 256>>>(Frac_verts_device_ptr,          // the pointer to device vector
                                                                          TYR,                            // seed
                                                                          DSIZE,                          // number of fracture
                                                                          L,                              // domain size
                                                                          ModeOfFractureSizeDistribution, // distribution pattern of fracture sizes
                                                                          SizeDistributionParameters,     // parameters of distribution of fracture sizes
                                                                          kappa,                          // kappa value of fisher distribution
                                                                          conductivity_beta,              // beta
                                                                          conductivity_gamma,             // gamma
                                                                          DomainDimensionRatio);          // ratio

                cudaDeviceSynchronize();
                Frac_verts_host = Frac_verts_device; // copy data from device to host
                std::map<pair<size_t, size_t>, pair<cuDFNsys::Vector3<_DataType_>, cuDFNsys::Vector3<_DataType_>>> Intersection_map;

                cuDFNsys::IdentifyIntersection<_DataType_> identifyInters{Frac_verts_host.size(),
                                                                          Frac_verts_device_ptr,
                                                                          false,
                                                                          Intersection_map};

                std::vector<std::vector<size_t>> ListClusters;
                std::vector<size_t> Percolation_cluster;
                cuDFNsys::Graph<_DataType_> G{(size_t)DSIZE, Intersection_map};
                G.UseDFS(ListClusters);
                cuDFNsys::IdentifyPercolationCluster<_DataType_> IdentiClu{ListClusters,
                                                                           Frac_verts_host,
                                                                           perco_dir,
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

                cuDFNsys::IdentifyIntersection<_DataType_> identifyInters2{Frac_verts_host.size(),
                                                                           Frac_verts_device_ptr, true,
                                                                           Intersection_map};

                cuDFNsys::Graph<_DataType_> G2{(size_t)DSIZE, Intersection_map};
                G2.UseDFS(ListClusters);
                cuDFNsys::IdentifyPercolationCluster<_DataType_> IdentiClu2{ListClusters,
                                                                            Frac_verts_host, perco_dir,
                                                                            Percolation_cluster};

                Frac_verts_device.clear();
                Frac_verts_device.shrink_to_fit();

                cout << "DFN generated. Running time: " << cuDFNsys::CPUSecond() - iStart_DFN << " sec\n";

                if (IfFlowSimulation && Percolation_cluster.size() > 0)
                {

                    double iStart_sortingFractures = cuDFNsys::CPUSecond();

                    std::vector<size_t> Fracs_percol; // will store ID / No of fractures in percolating cluster

                    cuDFNsys::GetAllPercolatingFractures GetPer{Percolation_cluster,
                                                                ListClusters,
                                                                Fracs_percol};
                    // this function is simple, just collecting fractures in the percolating clusters

                    std::vector<pair<int, int>> IntersectionPair_percol;
                    // will store the intersection pair of fractures in percolation clusters

                    cuDFNsys::RemoveDeadEndFrac<_DataType_> RDEF{Fracs_percol,            // fractures' ID in percolating cluster
                                                                 IntersectionPair_percol, // intersection pair
                                                                 (size_t)perco_dir,       // flow direction
                                                                 Frac_verts_host,         // host vector of fractures
                                                                 Intersection_map,        // map of intersection
                                                                 false};                  // does not remove dead-end fractures; true = remove dead ends
                    // the above function just sorts fractures and does not remove dead end fractures

                    cout << "Running time of sorting fractures: " << cuDFNsys::CPUSecond() - iStart_sortingFractures << " sec\n";

                    double iStart_mesh = cuDFNsys::CPUSecond();

                    cout << "Meshing ..." << endl;

                    cuDFNsys::Mesh<_DataType_> mesh{Frac_verts_host,         // host vector of fractures, after removing fractures of dead end
                                                    IntersectionPair_percol, // intersection pair
                                                    &Fracs_percol,           // fractures' ID in percolating cluster
                                                    minGrid,                 // minimum grid size
                                                    maxGrid,                 // maximum grid size
                                                    perco_dir,               // flow direction
                                                    L,                       // domain size
                                                    DomainDimensionRatio};

                    cout << "Mesh finished. Running time: " << cuDFNsys::CPUSecond() - iStart_mesh << " sec\n";

                    cout << "MHFEM ing ..." << endl;
                    double iStart_mhfem = cuDFNsys::CPUSecond();

                    cuDFNsys::MHFEM<_DataType_> fem{mesh,            // mesh object
                                                    Frac_verts_host, // fractures
                                                    L,               // hydraulic head of inlet = 100 m
                                                    0,               // hydraulic head of outlet = 20 m
                                                    perco_dir,       // flow direction
                                                    L,               // domain size
                                                    DomainDimensionRatio};

                    if (fem.QError > 1 || isnan(fem.Permeability) == 1)
                    {
                        string AS = "Found large error or isnan for permeability, the error: " + to_string(fem.QError) + ", the permeability: " + to_string(fem.Permeability);
                        throw cuDFNsys::ExceptionsIgnore(AS);
                    }

                    Permeability_A[0] = fem.Permeability;
                    Q_error_A[0] = fem.QError;

                    cout << "MHFEM finished. Running time: " << cuDFNsys::CPUSecond() - iStart_mhfem << " sec\n";
                    //---------------------
                }

                cout << "Outputing data...\n\n";
                string groupname = "group_" + cuDFNsys::ToStringWithWidth(i + 1, 3);
                string DirStr = "_x";

                if (perco_dir == 1)
                    DirStr = "_y";
                else if (perco_dir == 2)
                    DirStr = "_z";

                vector<string> datasetname = {
                    "P33_total",
                    "P33_connected" + DirStr,
                    "Ratio_of_P33" + DirStr,
                    "P33_largest_cluster",
                    "P32_total",
                    "P32_connected" + DirStr,
                    "Ratio_of_P32" + DirStr,
                    "P32_largest_cluster",
                    "P30",
                    "P30_connected" + DirStr,
                    "Ratio_of_P30" + DirStr,
                    "P30_largest_cluster",
                    "Percolation_probability" + DirStr,
                    "n_I",
                    "Permeability" + DirStr,
                    "Q_error" + DirStr};

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

                vector<uint2> dim_ss(data_input.size(), make_uint2(1, 1));

                h5out.AddDatasetsWithOneGroup(filename, groupname,
                                              datasetname, data_input, dim_ss);

                _DataType_ i_p[1] = {(_DataType_)i + 1};
                if (i == 0)
                    h5out.AddDataset(filename, "N", "Loop_times", i_p, dim_e);
                else
                    h5out.OverWrite(filename, "N", "Loop_times", i_p, dim_e);

                //---------try finished
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
        }
    }
    catch (...)
    {
        throw;
        return 0;
    };

    return 0;
};