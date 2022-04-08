///////////////////////////////////////////////////////////////////
// NAME:              RemoveDeadEndFrac.cuh
//
// PURPOSE:           Remove dead end fractures
//                    that do not conduct flow
//
// FUNCTIONS/OBJECTS: RemoveDeadEndFrac
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../GlobalDef/GlobalDef.cuh"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Fracture.cuh"
#include <vector>

using namespace std;
using namespace Eigen;

namespace cuDFNsys
{
class RemoveDeadEndFrac
{
public:
    RemoveDeadEndFrac(std::vector<size_t> &One_cluster,
                      std::vector<pair<int, int>> &Intersection_pair,
                      const size_t &dir,
                      const thrust::host_vector<cuDFNsys::Fracture> &Fracs,
                      std::map<pair<size_t, size_t>, pair<float3, float3>> Intersection_map);
};
}; // namespace cuDFNsys