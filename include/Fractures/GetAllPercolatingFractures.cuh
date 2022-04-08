///////////////////////////////////////////////////////////////////
// NAME:              GetAllPercolatingFractures.cuh
//
// PURPOSE:           Get all fractures in percolating clusters
//
// FUNCTIONS/OBJECTS: GetAllPercolatingFractures
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../GlobalDef/GlobalDef.cuh"
namespace cuDFNsys
{

class GetAllPercolatingFractures
{
public:
    GetAllPercolatingFractures(const std::vector<size_t> &Percolation_cluster,
                               const std::vector<std::vector<size_t>> &ListClusters,
                               std::vector<size_t> &Fracs_percol);
};

}; // namespace cuDFNsys