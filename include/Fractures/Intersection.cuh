///////////////////////////////////////////////////////////////////
// NAME:              Intersection.cuh
//
// PURPOSE:           Intersection struct
//
// FUNCTIONS/OBJECTS: Intersection
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../GlobalDef/GlobalDef.cuh"

namespace cuDFNsys
{
struct Intersection
{
    int2 FracIDPair;
    float3 Coord[2];
};
}; // namespace cuDFNsys