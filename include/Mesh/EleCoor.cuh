///////////////////////////////////////////////////////////////////
// NAME:              EleCoor.cuh
//
// PURPOSE:           Element coordinates
//
// FUNCTIONS/OBJECTS: EleCoor
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../GlobalDef/GlobalDef.cuh"

namespace cuDFNsys
{
struct EleCoor
{
    float x[3];
    float y[3];
};
}; // namespace cuDFNsys