// ====================================================
// NAME:        User's interface
// DESCRIPTION: Call cuDFNsys functions to do simulation.
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================

#include "cuDFNsys.cuh"
#include <unistd.h>
int main(int argc, char *argv[])
{
    int iDev = 0;
    GpuErrCheck(cudaSetDevice(iDev));
    int DSIZE = atoi(argv[1]);
    float L = atof(argv[2]);
    cuDFNsys::Warmup<<<DSIZE / 256 + 1, 256 /*  1, 2*/>>>();
    cudaDeviceSynchronize();
    cout << "Warmup finished\n";

    // cuDFNsys::Fracture sf;
    thrust::host_vector<cuDFNsys::Fracture> Frac_verts_host(DSIZE);
    thrust::device_vector<cuDFNsys::Fracture> Frac_verts_device(DSIZE);
    time_t t;
    time(&t);
    cuDFNsys::Fracture *Frac_verts_device_ptr = thrust::raw_pointer_cast(Frac_verts_device.data());

    Fractures<<<DSIZE / 256 + 1, 256 /*  1, 2*/>>>(Frac_verts_device_ptr,
                                                   (unsigned long)t,
                                                   DSIZE,
                                                   L,
                                                   1.5f,
                                                   1.0f,
                                                   15.0f,
                                                   0.0f,
                                                   0.0f);

    cudaDeviceSynchronize();
    Frac_verts_host = Frac_verts_device;
    cout << "DFN 1\n";
    double TStart = cuDFNsys::CpuSecond();
    std::map<pair<size_t, size_t>, pair<float3, float3>> Intersection_map;
    cuDFNsys::IdentifyIntersection identifyInters{Frac_verts_host, false, Intersection_map};
    double TElapse = cuDFNsys::CpuSecond() - TStart;
    cout << "running time of DFN 1 intersection: " << TElapse << "s \n";

    std::vector<std::vector<size_t>> ListClusters;
    std::vector<size_t> Percolation_cluster;
    TStart = cuDFNsys::CpuSecond();
    cuDFNsys::Graph G{(size_t)DSIZE, Intersection_map};
    G.UseDFS(ListClusters);
    TElapse = cuDFNsys::CpuSecond() - TStart;
    cout << "running time of DFN 1 DFS: " << TElapse << "s \n";

    cuDFNsys::IdentifyPercolationCluster IdentiClu{ListClusters, Frac_verts_host, 2, Percolation_cluster};
    cuDFNsys::MatlabPlotDFN As{"DFN_I.mat", "DFN_I.m", Frac_verts_host, Intersection_map, ListClusters, Percolation_cluster, false, true, true, true, L, 2};

    Intersection_map.clear();
    ListClusters.clear();
    Percolation_cluster.clear();
    cuDFNsys::IdentifyIntersection identifyInters2{Frac_verts_host, true, Intersection_map};

    cuDFNsys::Graph G2{(size_t)DSIZE, Intersection_map};
    G2.UseDFS(ListClusters);
    cuDFNsys::IdentifyPercolationCluster IdentiClu2{ListClusters, Frac_verts_host, 2, Percolation_cluster};
    cuDFNsys::MatlabPlotDFN Ak{"DFN_II.mat", "DFN_II.m", Frac_verts_host, Intersection_map, ListClusters, Percolation_cluster, true, true, true, true, L, 2};

    if (Percolation_cluster.size() > 0)
    {
        std::vector<size_t> Fracs_percol;
        cuDFNsys::GetAllPercolatingFractures GetPer{Percolation_cluster,
                                                    ListClusters,
                                                    Fracs_percol};
        std::vector<pair<int, int>> IntersectionPair_percol;
        cuDFNsys::RemoveDeadEndFrac RDEF{Fracs_percol,
                                         IntersectionPair_percol,
                                         2,
                                         Frac_verts_host,
                                         Intersection_map};
        cout << "\tmeshing ..." << endl;
        //cudaDeviceReset();
        cuDFNsys::Mesh meshr{Frac_verts_host, IntersectionPair_percol,
                             &Fracs_percol, 1.5, 2, 2, L};
        meshr.MatlabPlot("DFNMesh.mat",
                         "DFNMesh.m", Frac_verts_host, L, false, true);
    };
    return 0;
};