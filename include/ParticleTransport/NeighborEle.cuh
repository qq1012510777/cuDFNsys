///////////////////////////////////////////////////////////////////
// NAME:              NeighborEle.cuh
//
// PURPOSE:           The neighbor elements of an element
//
// FUNCTIONS/OBJECTS: NeighborEle
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../GlobalDef/GlobalDef.cuh"
#include "../DataTypeSelector/DataTypeSelector.cuh"

namespace cuDFNsys
{
struct NeighborEle
{
    // number of shared elements for an edge, (at least = 1, at most = _NumOfSharedEleAtMost in GlobalDef.cuh)
    uint NumNeighborEle = 0;
    // here, the default of number of shared elements is _NumOfNeighborEleAtMost.
    uint EleID[_NumOfNeighborEleAtMost] = {0};
};
}; // namespace cuDFNsys