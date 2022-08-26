///////////////////////////////////////////////////////////////////
// NAME:              Particle.cuh
//
// PURPOSE:           a struct of Particle: record the position of
//                    a particle
//
// FUNCTIONS/OBJECTS: Particle
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../DataTypeSelector/DataTypeSelector.cuh"
#include "../GlobalDef/GlobalDef.cuh"

namespace cuDFNsys
{
template <typename T>
struct Particle
{
    // position in 2d
    cuDFNsys::Vector2<T> Position2D;
    // element ID
    uint ElementID;
    // particle ID
    int ParticleID; // ParticleID == -1 means that this particle reaches outlet plane
};
}; // namespace cuDFNsys