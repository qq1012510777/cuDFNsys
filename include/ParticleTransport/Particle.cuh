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
#include "../MHFEM/MHFEM.cuh"

namespace cuDFNsys
{
struct Particle
{
    // position in 2d
    float2 Position2D;
    // element ID
    uint ElementID;
};
}; // namespace cuDFNsys