///////////////////////////////////////////////////////////////////
// NAME:              IdentifyPercolationCluster.cuh
//
// PURPOSE:           Identify Percolation Cluster,
//
// FUNCTIONS/OBJECTS: IdentifyPercolationCluster
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////

#pragma once
#include "../DataTypeSelector/DataTypeSelector.cuh"
#include "../GlobalDef/GlobalDef.cuh"
#include "Fracture.cuh"

namespace cuDFNsys
{
template <typename T>
class IdentifyPercolationCluster
{
public:
    // constructor
    IdentifyPercolationCluster(const std::vector<std::vector<size_t>> &ListClusters,
                               const thrust::host_vector<cuDFNsys::Fracture<T>> &Frac_verts_host,
                               const int &dir,
                               std::vector<size_t> &Percolation_cluster);
};
}; // namespace cuDFNsys