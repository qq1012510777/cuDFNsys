///////////////////////////////////////////////////////////////////
// NAME:              IdentifyIntersection.cuh
//
// PURPOSE:           Identify intersection
//
// FUNCTIONS/OBJECTS: IdentifyIntersection
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../Geometry/Geometry.cuh"
#include "../GlobalDef/GlobalDef.cuh"
#include "../MatrixManipulation/MatrixManipulation.cuh"
#include "Fracture.cuh"
#include "IdentifyFracPairSphericalDetection.cuh"
#include "IdentifyIntersectionKernel.cuh"
#include "Intersection.cuh"
#include <unistd.h>

namespace cuDFNsys
{
template <typename T>
class IdentifyIntersection
{
public:
    // constructor CPU
    IdentifyIntersection(thrust::host_vector<cuDFNsys::Fracture<T>> verts,
                         const bool &if_trucncated,
                         std::map<pair<size_t, size_t>, pair<cuDFNsys::Vector3<T>, cuDFNsys::Vector3<T>>> &Intersection_map);
    // constructor GPU
    IdentifyIntersection(const size_t &Fracsize,
                         cuDFNsys::Fracture<T> *Frac_verts_device_ptr,
                         const bool &if_trucncated,
                         std::map<pair<size_t, size_t>, pair<cuDFNsys::Vector3<T>, cuDFNsys::Vector3<T>>> &Intersection_map);
};
}; // namespace cuDFNsys