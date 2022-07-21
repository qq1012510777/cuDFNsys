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
#include "../DataTypeSelector/DataTypeSelector.cuh"
#include "../GlobalDef/GlobalDef.cuh"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Fracture.cuh"
#include <vector>

using namespace std;
using namespace Eigen;

namespace cuDFNsys
{
template <typename T>
class RemoveDeadEndFrac
{
public:
    RemoveDeadEndFrac(std::vector<size_t> &One_cluster,
                      std::vector<pair<int, int>> &Intersection_pair,
                      const size_t &dir,
                      thrust::host_vector<cuDFNsys::Fracture<T>> &Fracs,
                      std::map<pair<size_t, size_t>, pair<cuDFNsys::Vector3<T>, cuDFNsys::Vector3<T>>> Intersection_map);
};
}; // namespace cuDFNsys