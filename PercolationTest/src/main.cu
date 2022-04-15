// ====================================================
// NAME:        Percolation
// DESCRIPTION: My percolation program.
// AUTHOR:      Tingchang YIN
// DATE:        09/04/2022
// ====================================================

#include "cuDFNsys.cuh"
#include <unistd.h>
int main(int argc, char *argv[])
{
    int iDev = 0;
    GPUErrCheck(cudaSetDevice(iDev));

    double iStart = cuDFNsys::CPUSecond();

    uint ModeSizeDistr = (uint)(atoi(argv[1]));
    float4 SizeParaDis = make_float4(atof(argv[2]),
                                     atof(argv[3]),
                                     atof(argv[4]),
                                     atof(argv[5]));

    float kappa = atof(argv[6]);
    float conductivity_powerlaw_e = atof(argv[7]);

    float L = atof(argv[8]);
    int perco_dir = atoi(argv[9]);

    int Init_NUM_Frac = atoi(argv[10]);

    int Frac_increment = atoi(argv[11]);

    int LOOP_times = atoi(argv[12]);
    int inti_LOOP_times = 0;

    float minGrid = atof(argv[13]);
    float maxGrid = atof(argv[14]);
    float InletP = atof(argv[15]);
    float OutletP = atof(argv[16]);

    int MODEL_NO = atoi(argv[17]);
    int MC_NO = atoi(argv[18]);

    float P33_total_A[1] = {0};
    float P33_connected_A[1] = {0};
    float Ratio_of_P33_A[1] = {0};
    float P33_largest_cluster_A[1] = {0};

    float P32_total_A[1] = {0};
    float P32_connected_A[1] = {0};
    float Ratio_of_P32_A[1] = {0};
    float P32_largest_cluster_A[1] = {0};

    float P30_A[1] = {0};
    float P30_connected_A[1] = {0};
    float Ratio_of_P30_A[1] = {0};
    float P30_largest_cluster_A[1] = {0};

    float Percolation_probability_A[1] = {0};
    float n_I_A[1] = {0};

    float Permeability_A[1] = {0};
    float Q_error_A[1] = {0};

    cuDFNsys::HDF5API h5out;
    string filename = "Datafile_Model_" + cuDFNsys::ToStringWithWidth(MODEL_NO, 3) +
                      "_MC_" + cuDFNsys::ToStringWithWidth(MC_NO, 5) + ".h5";

    uint2 dim_e = make_uint2(1, 1);

    if (!h5out.IfH5FileExist(filename))
    {
        h5out.NewFile(filename);
        float domainsize[1] = {L};
        h5out.AddDataset(filename, "N", "Domain_size", domainsize, dim_e);
    }
    else
    {
        vector<double> Looptimes_k =
            h5out.ReadDataset(filename, "N", "Loop_times");

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
    for (int i = inti_LOOP_times; i < LOOP_times; ++i)
    {
    ReGen_111:;
        try
        {
            cout << "LOOP " << i << endl;
            double iStart_DFN = cuDFNsys::CPUSecond();

            int DSIZE = Init_NUM_Frac + (i + 1) * Frac_increment;

            thrust::host_vector<cuDFNsys::Fracture> Frac_verts_host(DSIZE);
            thrust::device_vector<cuDFNsys::Fracture> Frac_verts_device(DSIZE);
            cuDFNsys::Fracture *Frac_verts_device_ptr = thrust::raw_pointer_cast(Frac_verts_device.data());
            //cout << "\tAllocated memory to vectors of cuDFNsys::Fracture\n";
            cuDFNsys::Warmup<<<DSIZE / 256 + 1, 256 /*  1, 2*/>>>();
            cudaDeviceSynchronize();
            time_t t;
            time(&t);

            cuDFNsys::Fractures<<<DSIZE / 256 + 1, 256 /*  1, 2*/>>>(Frac_verts_device_ptr,
                                                                     (unsigned long)t,
                                                                     DSIZE,
                                                                     L,
                                                                     ModeSizeDistr,
                                                                     SizeParaDis,
                                                                     conductivity_powerlaw_e);
            //-----------
            cudaDeviceSynchronize();
            //cout << "\tGenerated fractures\n";
            Frac_verts_host = Frac_verts_device;

            std::map<pair<size_t, size_t>, pair<float3, float3>> Intersection_map;
            try
            {
                cuDFNsys::IdentifyIntersection identifyInters{Frac_verts_host.size(),
                                                              Frac_verts_device_ptr,
                                                              false,
                                                              Intersection_map};
            }
            catch (thrust::system::system_error &e)
            {
                cout << "\033[1;32m\tThrew an exception: " << e.what() << "\n\tLet's use CPU finish identification of intersections\033[0m" << endl;
                cudaDeviceReset();
                Intersection_map.clear();
                cuDFNsys::IdentifyIntersection identifyInters3{Frac_verts_host,
                                                               false,
                                                               Intersection_map};
            };
            //cout << "\tIdentified intersection of DFN with complete factures\n";
            std::vector<std::vector<size_t>> ListClusters;
            std::vector<size_t> Percolation_cluster;
            cuDFNsys::Graph G{(size_t)DSIZE, Intersection_map};
            G.UseDFS(ListClusters);
            cuDFNsys::IdentifyPercolationCluster IdentiClu{ListClusters,
                                                           Frac_verts_host, perco_dir,
                                                           Percolation_cluster};

            Intersection_map.clear();
            ListClusters.clear();
            Percolation_cluster.clear();
            try
            {
                cuDFNsys::IdentifyIntersection identifyInters2{Frac_verts_host.size(),
                                                               Frac_verts_device_ptr, true,
                                                               Intersection_map};
            }
            catch (thrust::system::system_error &e)
            {
                cout << "\033[1;32m\tThrew an exception: " << e.what() << "\n\tLet's use CPU finish identification of intersections\033[0m" << endl;
                cudaDeviceReset();
                Intersection_map.clear();
                cuDFNsys::IdentifyIntersection identifyInters3{Frac_verts_host,
                                                               true,
                                                               Intersection_map};
            };
            //cout << "\tIdentified intersection of DFN with truncated factures\n";
            cuDFNsys::Graph G2{(size_t)DSIZE, Intersection_map};
            G2.UseDFS(ListClusters);
            cuDFNsys::IdentifyPercolationCluster IdentiClu2{ListClusters,
                                                            Frac_verts_host, perco_dir,
                                                            Percolation_cluster};
            cuDFNsys::GetStatistics(Frac_verts_host,
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

            double iElaps_DFN = cuDFNsys::CPUSecond() - iStart_DFN;
            cout << "\tDFN generated. Running time: " << iElaps_DFN << " sec\n";

            Frac_verts_device.clear();
            Frac_verts_device.shrink_to_fit();
            //cudaDeviceReset();

            //Connections.clear();
            //Connections.shrink_to_fit();

            //-------------------------
            if (Percolation_cluster.size() > 0)
            {
                double iStart_mesh = cuDFNsys::CPUSecond();
                std::vector<size_t> Fracs_percol;
                cuDFNsys::GetAllPercolatingFractures GetPer{Percolation_cluster,
                                                            ListClusters,
                                                            Fracs_percol};

                std::vector<pair<int, int>> IntersectionPair_percol;

                cuDFNsys::RemoveDeadEndFrac RDEF{Fracs_percol,
                                                 IntersectionPair_percol,
                                                 (size_t)perco_dir,
                                                 Frac_verts_host,
                                                 Intersection_map};

                cout << "\tMeshing ..." << endl;
                //cudaDeviceReset();
                cuDFNsys::Mesh mesh{Frac_verts_host, IntersectionPair_percol,
                                    &Fracs_percol, minGrid, maxGrid, perco_dir, L};

                double iElaps_mesh = cuDFNsys::CPUSecond() - iStart_mesh;
                cout << "\tMesh finished. Running time: " << iElaps_mesh << " sec\n";

                cout << "\tMHFEM ing ..." << endl;
                double iStart_mhfem = cuDFNsys::CPUSecond();
                cuDFNsys::MHFEM fem{mesh, Frac_verts_host, InletP, OutletP, perco_dir, L};

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
            cudaDeviceReset();
            cout << "\tOutputing data...\n";
            string groupname = "group_" + cuDFNsys::ToStringWithWidth(i + 1, 3);
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
            vector<float *> data_input = {P33_total_A,
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
            float i_p[1] = {(float)i + 1};
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