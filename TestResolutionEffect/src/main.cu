// ====================================================
// NAME:        TestResolutionEffect
// DESCRIPTION:
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================

#include "cuDFNsys.cuh"
#include <unistd.h>
int main(int argc, char *argv[])
{

    try
    {
        double istart = cuDFNsys::CPUSecond();

        int dev = 0;
        GPUErrCheck(cudaSetDevice(dev));

        int DSIZE = 0;
        float L = 0;
        float minGrid = 0;
        float maxGrid = 0;

        DSIZE = atoi(argv[1]);
        L = atof(argv[2]);

        int perco_dir = 2;

        cout << "preparing" << endl;

        thrust::host_vector<cuDFNsys::Fracture> Frac_verts_host(DSIZE);
        thrust::device_vector<cuDFNsys::Fracture> Frac_verts_device(DSIZE);
        cuDFNsys::Fracture *Frac_verts_device_ptr;
        Frac_verts_device_ptr = thrust::raw_pointer_cast(Frac_verts_device.data());

        cuDFNsys::Warmup<<<DSIZE / 256 + 1, 256 /*  1, 2*/>>>();
        cudaDeviceSynchronize();
        time_t t;
        time(&t);

        cout << "generating fractures" << endl;

        cuDFNsys::Fractures<<<DSIZE / 256 + 1, 256 /*  1, 2*/>>>(Frac_verts_device_ptr,
                                                                 (unsigned long)t,
                                                                 DSIZE,
                                                                 L,
                                                                 atof(argv[3]),
                                                                 atof(argv[4]),
                                                                 atof(argv[5]),
                                                                 atof(argv[6]),
                                                                 atof(argv[7]));
        cudaDeviceSynchronize();
        Frac_verts_host = Frac_verts_device;
        cout << "identifying intersections with complete fractures" << endl;
        std::map<pair<size_t, size_t>, pair<float3, float3>> Intersection_map;
        cuDFNsys::IdentifyIntersection identifyInters{Frac_verts_host.size(),
                                                      Frac_verts_device_ptr,
                                                      false,
                                                      Intersection_map};
        cout << "identifying cluster with complete fractures" << endl;
        std::vector<std::vector<size_t>> ListClusters;
        std::vector<size_t> Percolation_cluster;
        cuDFNsys::Graph G{(size_t)DSIZE, Intersection_map};
        G.UseDFS(ListClusters);
        cuDFNsys::IdentifyPercolationCluster IdentiClu{ListClusters,
                                                       Frac_verts_host, perco_dir,
                                                       Percolation_cluster};
        cout << "DFN I finished" << endl;
        cuDFNsys::MatlabPlotDFN As{"DFN_I.mat", "DFN_I.m",
                                   Frac_verts_host, Intersection_map, ListClusters,
                                   Percolation_cluster, false, true, true, true,
                                   L, perco_dir};
        //
        Intersection_map.clear();
        ListClusters.clear();
        Percolation_cluster.clear();
        cout << "identifying intersections with truncated fractures" << endl;
        cuDFNsys::IdentifyIntersection identifyInters2{Frac_verts_host.size(),
                                                       Frac_verts_device_ptr, true,
                                                       Intersection_map};
        cout << "identifying cluster with truncated fractures" << endl;
        cuDFNsys::Graph G2{(size_t)DSIZE, Intersection_map};
        G2.UseDFS(ListClusters);
        cuDFNsys::IdentifyPercolationCluster IdentiClu2{ListClusters,
                                                        Frac_verts_host, perco_dir,
                                                        Percolation_cluster};
        cout << "DFN II finished" << endl;
        cuDFNsys::MatlabPlotDFN As2{"DFN_II.mat", "DFN_II.m",
                                    Frac_verts_host, Intersection_map, ListClusters,
                                    Percolation_cluster, true, true, true, true,
                                    L, perco_dir};
        Frac_verts_device.clear();
        Frac_verts_device.shrink_to_fit();
        //-----------
        if (Percolation_cluster.size() > 0)
        {
            for (size_t i = 0;; ++i)
            {
                cout << "Please input minimum and maximum grid sizes:\n";
                cin >> minGrid;
                cin >> maxGrid;
                cout << "minGrid: " << minGrid << ", " << maxGrid << endl;
                double istart_1 = cuDFNsys::CPUSecond();
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
                cout << "meshing ..." << endl;

                cuDFNsys::Mesh mesh{Frac_verts_host, IntersectionPair_percol,
                                    &Fracs_percol, minGrid, maxGrid, perco_dir, L};

                cout << "MHFEM ing ..." << endl;

                cuDFNsys::MHFEM fem{mesh, Frac_verts_host, 100, 20, perco_dir, L};
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

                mesh.MatlabPlot("DFN_mesh_" + to_string(i + 1) + ".mat",
                                "DFN_mesh_" + to_string(i + 1) + ".m",
                                Frac_verts_host, L, true, true);
                fem.MatlabPlot("MHFEM_" + to_string(i + 1) + ".mat",
                               "MHFEM_" + to_string(i + 1) + ".m",
                               mesh, L);

                int uy = 0;
                cout << "keep changing grid size?\n";
                cin >> uy;
                if (uy == 0)
                    break;
            };
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