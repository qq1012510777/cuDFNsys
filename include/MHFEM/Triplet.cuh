///////////////////////////////////////////////////////////////////
// NAME:              Triplet.cuh
//
// PURPOSE:           triplet to generate sparse matrix
//
// FUNCTIONS/OBJECTS: Triplet
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../GlobalDef/GlobalDef.cuh"
#include "../DataTypeSelector/DataTypeSelector.cuh"

namespace cuDFNsys
{
// triplet to generate sparse matrix
template <typename T>
struct Triplet
{
    int row = -1;
    int col = -1;
    T val = -1;
};
}; // namespace cuDFNsys
