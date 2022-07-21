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
#include "../DataTypeSelector/DataTypeSelector.cuh"
#include "../GlobalDef/GlobalDef.cuh"

namespace cuDFNsys
{
template <typename T>
struct EleCoor
{
    T x[3];
    T y[3];
};
}; // namespace cuDFNsys