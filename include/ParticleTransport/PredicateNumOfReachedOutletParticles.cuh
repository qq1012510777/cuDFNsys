///////////////////////////////////////////////////////////////////
// NAME:              PredicateNumOfReachedOutletParticles.cuh
//
// PURPOSE:           a predicate for thrust::count_if to 
//                    count how many particles reach outlet plane
//
// FUNCTIONS/OBJECTS: PredicateNumOfReachedOutletParticles
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../DataTypeSelector/DataTypeSelector.cuh"
#include "../GlobalDef/GlobalDef.cuh"
#include "Particle.cuh"

namespace cuDFNsys
{
template <typename T>
struct PredicateNumOfReachedOutletParticles
{
    __host__ __device__ bool operator()(const cuDFNsys::Particle<T> &x) const
    {
        return (x.ParticleID == -1 ? true : false);
    };
};
}; // namespace cuDFNsys