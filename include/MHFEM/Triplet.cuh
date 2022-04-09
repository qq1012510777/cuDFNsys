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

namespace cuDFNsys
{
// triplet to generate sparse matrix
struct Triplet
{
    int row = -1;
    int col = -1;
    float val = -1;
};
}; // namespace cuDFNsys
