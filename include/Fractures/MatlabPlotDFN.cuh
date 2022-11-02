///////////////////////////////////////////////////////////////////
// NAME:              MatlabPlotDFN.cuh
//
// PURPOSE:           Plot DFN with m and mat file
//                    m is command file; mat is data file
//                    updated: this function can also output py and h5 file
//
// FUNCTIONS/OBJECTS: MatlabPlotDFN
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../DataTypeSelector/DataTypeSelector.cuh"
#include "../GlobalDef/GlobalDef.cuh"
#include "../HDF5API/HDF5API.cuh"
//#include "../MatlabAPI/MatlabAPI.cuh"
#include "Fracture.cuh"
#include <algorithm>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace cuDFNsys
{
template <typename T>
class MatlabPlotDFN
{
public:
    // constructor
    MatlabPlotDFN(string mat_key,                                                                                     // mat file name
                  string command_key,                                                                                 // m file name
                  thrust::host_vector<cuDFNsys::Fracture<T>> Frac_verts_host,                                         // Vector of Fracture
                  std::map<pair<size_t, size_t>, pair<cuDFNsys::Vector3<T>, cuDFNsys::Vector3<T>>> Intersection_host, // vector of Intersection
                  std::vector<std::vector<size_t>> ListClusters,                                                      // List of clusters: a cluster is a vector
                  std::vector<size_t> Percolation_cluster,                                                            // Percolation cluster contains the cluster vector ID
                  bool If_show_truncated_frac,                                                                        // if show truncated fractures
                  bool If_show_intersection,
                  bool If_show_cluster,
                  bool If_show_orientation,
                  T L,
                  int dir,
                  bool if_python_visualization = false,
                  string PythonName_Without_suffix = "DFN_py");
};
}; // namespace cuDFNsys