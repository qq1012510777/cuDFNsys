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
#include "../GlobalDef/GlobalDef.cuh"
#include "../DataTypeSelector/DataTypeSelector.cuh"

namespace cuDFNsys
{
template <typename T>
struct Particle
{
    // position in 2d
    cuDFNsys::Vector2<T> Position2D;
    // element ID
    uint ElementID;

    bool IfReachOutletPlane = false;
};
}; // namespace cuDFNsys