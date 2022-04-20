///////////////////////////////////////////////////////////////////
// NAME:              EdgeToEle.cuh
//
// PURPOSE:           a struct of EdgeToEle: record the shared element IDs of each
//                    (separated) edge
//
// FUNCTIONS/OBJECTS: EdgeToEle
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../MHFEM/MHFEM.cuh"

namespace cuDFNsys
{
struct EdgeToEle
{
    // number of shared elements for an edge, (at least = 1, at most = _NumOfSharedEleAtMost in GlobalDef.cuh)
    uint NumSharedEle = 0;
    // here, the default of number of shared elements is _NumOfSharedEleAtMost.
    uint EleID[_NumOfSharedEleAtMost] = {0};
};
}; // namespace cuDFNsys