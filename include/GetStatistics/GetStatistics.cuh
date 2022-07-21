///////////////////////////////////////////////////////////////////
// NAME:              GetStatistics.h
//
// PURPOSE:           Get P30, P32, permeability and so on
//
// FUNCTIONS/OBJECTS: GetStatistics
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../DataTypeSelector/DataTypeSelector.cuh"
#include "../Fractures/Fracture.cuh"
#include "../GlobalDef/GlobalDef.cuh"

namespace cuDFNsys
{
template <typename T>
void GetStatistics(const thrust::host_vector<cuDFNsys::Fracture<T>> &Frac_verts_host,
                   const std::map<pair<size_t, size_t>, pair<cuDFNsys::Vector3<T>, cuDFNsys::Vector3<T>>> &Intersection_map,
                   const std::vector<std::vector<size_t>> &ListClusters,
                   const std::vector<size_t> &Percolation_cluster,
                   const T L,
                   T &P33_total_B,
                   T &P33_connected_B,
                   T &Ratio_of_P33_B,
                   T &P33_largest_cluster_B,
                   T &P32_total_B,
                   T &P32_connected_B,
                   T &Ratio_of_P32_B,
                   T &P32_largest_cluster_B,
                   T &P30_B,
                   T &P30_connected_B,
                   T &Ratio_of_P30_B,
                   T &P30_largest_cluster_B,
                   T &Percolation_probability_B,
                   T &n_I_B);
}; // namespace cuDFNsys