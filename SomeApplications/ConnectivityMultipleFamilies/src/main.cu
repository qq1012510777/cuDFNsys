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
// DESCRIPTION: Generate multiple families and percolation test
// AUTHOR:      Tingchang YIN
// DATE:        04/03/2023
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

        uint NumFamilies = atof(argv[2]);

        uint count_Argv_index = 2;

        std::vector<cuDFNsys::Vector4<_DataType_>> ParaSizeDistri(NumFamilies);
        std::vector<_DataType_> kappa_(NumFamilies);
        std::vector<_DataType_> beta_(NumFamilies);
        std::vector<_DataType_> gamma_(NumFamilies);
        std::vector<cuDFNsys::Vector3<_DataType_>> rotationAxis(NumFamilies);
        std::vector<_DataType_> AngleRotation(NumFamilies);
        std::vector<uint> DSIZE_incre(NumFamilies);

        for (uint i = 0; i < NumFamilies; ++i)
        {
            char *FF = argv[count_Argv_index + 1];
            if (*FF != 'F' && *FF != 'f')
                throw cuDFNsys::ExceptionsPause("Input of the attributes of one fracture family should be started with a character 'f' or 'F'\n");

            count_Argv_index++;

            DSIZE_incre[i] = atoi(argv[count_Argv_index + 1]);
            count_Argv_index++;

            ParaSizeDistri[i] =
                cuDFNsys::MakeVector4((_DataType_)atof(argv[count_Argv_index + 1]),
                                      (_DataType_)atof(argv[count_Argv_index + 2]),
                                      (_DataType_)atof(argv[count_Argv_index + 3]),
                                      (_DataType_)atof(argv[count_Argv_index + 4]));

            count_Argv_index += 4;

            kappa_[i] = (_DataType_)atof(argv[count_Argv_index + 1]);
            beta_[i] = (_DataType_)atof(argv[count_Argv_index + 2]);
            gamma_[i] = (_DataType_)atof(argv[count_Argv_index + 3]);
            count_Argv_index += 3;

            rotationAxis[i] = cuDFNsys::MakeVector3(
                (_DataType_)atof(argv[count_Argv_index + 1]),
                (_DataType_)atof(argv[count_Argv_index + 2]),
                (_DataType_)atof(argv[count_Argv_index + 3]));
            count_Argv_index += 3;

            AngleRotation[i] =
                (_DataType_)atof(argv[count_Argv_index + 1]); // degree
            count_Argv_index++;
        }

        char *FF = argv[count_Argv_index + 1];

        if (*FF != 'D' && *FF != 'd')
            throw cuDFNsys::ExceptionsPause("Finishing the input of fracture attributes should end up with a character 'D' or 'd'\n");

        count_Argv_index++;

        uint LoopTimes = atoi(argv[count_Argv_index + 1]);
        count_Argv_index++;

        int MODEL_NO = atoi(argv[count_Argv_index + 1]);
        count_Argv_index++;

        int MC_NO = atoi(argv[count_Argv_index + 1]);
        count_Argv_index++;

        int inti_LOOP_times = 1;

        cout << "LoopTimes: " << LoopTimes << endl;
        cout << "NumFamilies: " << NumFamilies << endl;
        cout << "MODEL_NO: " << MODEL_NO << endl;
        cout << "MC_NO: " << MC_NO << endl;
        for (uint i = 0; i < NumFamilies; ++i)
        {
            cout << "Family NO. " << i + 1 << endl;
            cout << "\tParaSizeDistri: " << ParaSizeDistri[i].x << ", "
                 << ParaSizeDistri[i].y << ", "
                 << ParaSizeDistri[i].z << ", "
                 << ParaSizeDistri[i].w << endl;
            cout << "\tkappa_: " << kappa_[i] << endl;
            cout << "\tbeta_: " << beta_[i] << endl;
            cout << "\tgamma_: " << gamma_[i] << endl;
            cout << "\trotationAxis: " << rotationAxis[i].x << ", "
                 << rotationAxis[i].y << ", "
                 << rotationAxis[i].z << ", " << endl;
            cout << "\tAngleRotation: " << AngleRotation[i] << endl;
            cout << "\tDSIZE_incre: " << DSIZE_incre[i] << endl;
        }

        //--------------------------

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

                if (Looptimes_k[0] >= LoopTimes)
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
                inti_LOOP_times = 1;
            }
        }

        for (uint NP = inti_LOOP_times; NP <= LoopTimes; ++NP)
        {
        ReGen_111:;
            try
            {

                cout << "Loop: " << NP << endl;
                double iStart_DFN = cuDFNsys::CPUSecond();

                thrust::host_vector<cuDFNsys::Fracture<_DataType_>> Frac_verts_host;

                for (uint i = 0; i < NumFamilies; ++i)
                {
                    uint DSIZE = DSIZE_incre[i] * NP;

                    thrust::host_vector<cuDFNsys::Fracture<_DataType_>> Frac_verts_host_i(DSIZE);
                    thrust::device_vector<cuDFNsys::Fracture<_DataType_>> Frac_verts_device_i(DSIZE);
                    cuDFNsys::Fracture<_DataType_> *Frac_verts_device_ptr_i;
                    Frac_verts_device_ptr_i = thrust::raw_pointer_cast(Frac_verts_device_i.data());

                    srand((unsigned int)time(0) + i * 100);
                    Eigen::MatrixXd Ter = Eigen::MatrixXd::Random(1, 1);

                    // cout << "generating fracture group No. " << i + 1 << ", random seed: " << (unsigned long)ceil(abs(Ter(0, 0)) * ((unsigned long)t * 1.0)) << endl;
                    // cout << "\tnum of fractures: " << DSIZE << endl;
                    // cout << "\tangleRotation: " << AngleRotation * M_PI / 180.0 << endl;
                    // cout << "\trotationAxis: " << rotationAxis.x << ", " << rotationAxis.y << ", " << rotationAxis.z << endl;
                    cuDFNsys::Fractures<_DataType_><<<DSIZE / 256 + 1, 256>>>(Frac_verts_device_ptr_i,
                                                                              (unsigned long)t + (unsigned long)ceil(abs(Ter(0, 0)) * ((unsigned long)t * 1.0)),
                                                                              DSIZE, L,
                                                                              0, // power law
                                                                              ParaSizeDistri[i],
                                                                              kappa_[i],  // kappa
                                                                              beta_[i],   // beta
                                                                              gamma_[i]); // gamma
                    cudaDeviceSynchronize();
                    Frac_verts_host_i = Frac_verts_device_i;

                    //-------------rotate the fractures
                    if (AngleRotation[i] == 0)
                    {
                        Frac_verts_host.insert(Frac_verts_host.end(), Frac_verts_host_i.begin(), Frac_verts_host_i.end());
                        continue;
                    }

                    cuDFNsys::Quaternion<_DataType_> qua;
                    qua = qua.DescribeRotation(rotationAxis[i], AngleRotation[i] * M_PI / 180.0);

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

                uint DSIZE = Frac_verts_host.size(); // num of fractures
                cout << "\tNumFracs = " << DSIZE << endl;
                thrust::device_vector<cuDFNsys::Fracture<_DataType_>> Frac_verts_device = Frac_verts_host;
                cuDFNsys::Fracture<_DataType_> *Frac_verts_device_ptr;
                Frac_verts_device_ptr = thrust::raw_pointer_cast(Frac_verts_device.data());

                //cout << "dealing with all " << DSIZE << " fractures ..." << endl;
                cuDFNsys::FracturesChangeDomainSize<<<DSIZE / 256 + 1, 256>>>(Frac_verts_device_ptr,
                                                                              DSIZE, L);
                cudaDeviceSynchronize();
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

                double iElaps_DFN = cuDFNsys::CPUSecond() - iStart_DFN;
                cout << "\tDFN generated. Running time: " << iElaps_DFN << " sec\n";

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

                //-----------

                cout << "\tOutputing data...\n";
                string groupname = "group_" + cuDFNsys::ToStringWithWidth(NP, 3);

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
                    "n_I"};
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
                                                   n_I_A};
                uint2 dim_ki;
                dim_ki.x = 1;
                dim_ki.y = 1;
                vector<uint2> dim_ss(data_input.size(), dim_ki);
                h5out.AddDatasetsWithOneGroup(filename, groupname,
                                              datasetname, data_input, dim_ss);
                //cout << 1 << endl;
                _DataType_ i_p[1] = {(_DataType_)NP};
                if (NP == 1)
                    h5out.AddDataset(filename, "N", "Loop_times", i_p, dim_ki);
                else
                    h5out.OverWrite(filename, "N", "Loop_times", i_p, dim_ki);
            }
            catch (cuDFNsys::ExceptionsIgnore &e)
            {
                cout << e.what() << endl;
                goto ReGen_111;
            }
            catch (cuDFNsys::ExceptionsPause &e)
            {
                cout << e.what() << endl;
            }
            catch (...)
            {
                throw;
            };
        }

        cout << "Loop finished!\n";
        double iElaps = cuDFNsys::CPUSecond() - istart;
        cout << "Running times of this loop is: " << iElaps << " sec\n";
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