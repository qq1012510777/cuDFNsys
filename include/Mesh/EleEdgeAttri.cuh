///////////////////////////////////////////////////////////////////
// NAME:              EleEdgeAttri.cuh
//
// PURPOSE:           a struct of EleEdgeAttri
//                    edge NO and attributes
//
// FUNCTIONS/OBJECTS: EleEdgeAttri
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../GlobalDef/GlobalDef.cuh"

namespace cuDFNsys
{
struct EleEdgeAttri
{
    int e[3];  // 0: inlet, 1: outlet, 2: neumann, 3: interior
    int no[3]; // sep, sep, sep, global_interior
};
}; // namespace cuDFNsys