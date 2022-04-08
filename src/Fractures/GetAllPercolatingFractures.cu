#include "Fractures/GetAllPercolatingFractures.cuh"

// ====================================================
// NAME:        GetAllPercolatingFractures
// DESCRIPTION: Get all fractures in percolating clusters
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================
cuDFNsys::GetAllPercolatingFractures::GetAllPercolatingFractures(const std::vector<size_t> &Percolation_cluster,
                                                                 const std::vector<std::vector<size_t>> &ListClusters,
                                                                 std::vector<size_t> &Fracs_percol)
{
    size_t fracNUM = 0;

    for (size_t i = 0; i < Percolation_cluster.size(); ++i)
        fracNUM += ListClusters[Percolation_cluster[i]].size();

    Fracs_percol.clear();
    Fracs_percol.reserve(fracNUM);

    for (size_t i = 0; i < Percolation_cluster.size(); ++i)
        Fracs_percol.insert(Fracs_percol.end(),
                            ListClusters[Percolation_cluster[i]].begin(),
                            ListClusters[Percolation_cluster[i]].end());
}; // GetAllPercolatingFractures
