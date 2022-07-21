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
#include "../DataTypeSelector/DataTypeSelector.cuh"
#include "../GlobalDef/GlobalDef.cuh"

namespace cuDFNsys
{
template <typename T>
struct Intersection
{
    int2 FracIDPair;
    cuDFNsys::Vector3<T> Coord[2];
};
}; // namespace cuDFNsys