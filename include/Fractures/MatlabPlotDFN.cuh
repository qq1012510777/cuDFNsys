///////////////////////////////////////////////////////////////////
// NAME:              MatlabPlotDFN.cuh
//
// PURPOSE:           Plot DFN with m and mat file
//                    m is command file; mat is data file
//
// FUNCTIONS/OBJECTS: MatlabPlotDFN
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../GlobalDef/GlobalDef.cuh"
#include "../MatlabAPI/MatlabAPI.cuh"
#include "Fracture.cuh"
#include <algorithm>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace cuDFNsys
{
class MatlabPlotDFN
{
public:
    // constructor
    MatlabPlotDFN(string mat_key,                                          // mat file name
                  string command_key,                                      // m file name
                  thrust::host_vector<cuDFNsys::Fracture> Frac_verts_host, // Vector of Fracture
                  MapIntersection Intersection_host,                       // vector of Intersection
                  std::vector<std::vector<size_t>> ListClusters,           // List of clusters: a cluster is a vector
                  std::vector<size_t> Percolation_cluster,                 // Percolation cluster contains the cluster vector ID
                  bool If_show_truncated_frac,                             // if show truncated fractures
                  bool If_show_intersection,
                  bool If_show_cluster,
                  bool If_show_orientation,
                  float L,
                  int dir);
};
}; // namespace cuDFNsys