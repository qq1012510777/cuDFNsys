#include "Fractures/IdentifyPercolationCluster.cuh"

// ====================================================
// NAME:        IdentifyPercolationCluster
// DESCRIPTION: Identify Percolation Clusters
//              in a  DFN
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================
template <typename T>
cuDFNsys::IdentifyPercolationCluster<T>::IdentifyPercolationCluster(const std::vector<std::vector<size_t>> &ListClusters,
                                                                    const thrust::host_vector<cuDFNsys::Fracture<T>> &Frac_verts_host,
                                                                    const int &dir,
                                                                    std::vector<size_t> &Percolation_cluster)
{
    Percolation_cluster.clear();

    int inlet = 0, outlet = 0;
    if (dir == 0)
    {
        outlet = 1;
    }
    else if (dir == 1)
    {
        inlet = 2;
        outlet = 3;
    }
    else if (dir > 1)
    {
        inlet = 4;
        outlet = 5;
    }

    for (size_t i = 0; i < ListClusters.size(); ++i)
    {
        bool inlet_1 = false;
        bool outlet_1 = false;

        for (size_t j = 0; j < ListClusters[i].size(); ++j)
        {
            int FracID = ListClusters[i][j];

            if (Frac_verts_host[FracID].ConnectModelSurf[inlet] == true)
                inlet_1 = true;

            if (Frac_verts_host[FracID].ConnectModelSurf[outlet] == true)
                outlet_1 = true;

            if (inlet_1 == true && outlet_1 == true)
            {
                Percolation_cluster.push_back(i);
                break;
            }
        }
    };
    Percolation_cluster.shrink_to_fit();
}; // IdentifyPercolationCluster
template cuDFNsys::IdentifyPercolationCluster<double>::IdentifyPercolationCluster(const std::vector<std::vector<size_t>> &ListClusters,
                                                                                  const thrust::host_vector<cuDFNsys::Fracture<double>> &Frac_verts_host,
                                                                                  const int &dir,
                                                                                  std::vector<size_t> &Percolation_cluster);
template cuDFNsys::IdentifyPercolationCluster<float>::IdentifyPercolationCluster(const std::vector<std::vector<size_t>> &ListClusters,
                                                                                 const thrust::host_vector<cuDFNsys::Fracture<float>> &Frac_verts_host,
                                                                                 const int &dir,
                                                                                 std::vector<size_t> &Percolation_cluster);