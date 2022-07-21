///////////////////////////////////////////////////////////////////
// NAME:              DataTypeSelector.cuh
//
// PURPOSE:           define data type: double or float
//
// FUNCTIONS/OBJECTS: N/A
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////

#pragma once
#include "../GlobalDef/GlobalDef.cuh"

namespace cuDFNsys
{
template <typename T>
struct DataTypeSelector
{
};

template <>
struct DataTypeSelector<double>
{
    using Vector4Type = double4;
    using Vector3Type = double3;
    using Vector2Type = double2;
    using Vector1Type = double;
};

template <>
struct DataTypeSelector<float>
{
    using Vector4Type = float4;
    using Vector3Type = float3;
    using Vector2Type = float2;
    using Vector1Type = float;
};

template <typename P>
using Vector4 = typename DataTypeSelector<P>::Vector4Type;

template <typename P>
using Vector3 = typename DataTypeSelector<P>::Vector3Type;

template <typename P>
using Vector2 = typename DataTypeSelector<P>::Vector2Type;

template <typename P>
using Vector1 = typename DataTypeSelector<P>::Vector1Type;

}; // namespace cuDFNsys