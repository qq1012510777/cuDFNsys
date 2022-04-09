///////////////////////////////////////////////////////////////////
// NAME:              ToStringWithWidth.cuh
//
// PURPOSE:           convert number to string with fixed length
//
// FUNCTIONS/OBJECTS: ToStringWithWidth
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../GlobalDef/GlobalDef.cuh"

namespace cuDFNsys
{
template <class T>
string ToStringWithWidth(const T &val, const size_t &width);
}; // namespace cuDFNsys