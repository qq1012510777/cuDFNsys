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
#include "../Fractures/Fracture.cuh"
#include "../GlobalDef/GlobalDef.cuh"

namespace cuDFNsys
{
void GetStatistics(const thrust::host_vector<cuDFNsys::Fracture> &Frac_verts_host,
                   const std::map<pair<size_t, size_t>, pair<float3, float3>> &Intersection_map,
                   const std::vector<std::vector<size_t>> &ListClusters,
                   const std::vector<size_t> &Percolation_cluster,
                   const float L,
                   float &P33_total_B,
                   float &P33_connected_B,
                   float &Ratio_of_P33_B,
                   float &P33_largest_cluster_B,
                   float &P32_total_B,
                   float &P32_connected_B,
                   float &Ratio_of_P32_B,
                   float &P32_largest_cluster_B,
                   float &P30_B,
                   float &P30_connected_B,
                   float &Ratio_of_P30_B,
                   float &P30_largest_cluster_B,
                   float &Percolation_probability_B,
                   float &n_I_B);
}; // namespace cuDFNsys