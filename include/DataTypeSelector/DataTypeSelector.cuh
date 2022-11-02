/****************************************************************************
* cuDFNsys - simulating flow and transport in 3D fracture networks          *
* Copyright (C) 2022, Tingchang YIN, Sergio GALINDO-TORRES                  *
*                                                                           *
* This program is free software: you can redistribute it and/or modify      *
* it under the terms of the GNU Affero General Public License as            *
* published by the Free Software Foundation, either version 3 of the        *
* License, or (at your option) any later version.                           *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU Affero General Public License for more details.                       *
*                                                                           *
* You should have received a copy of the GNU Affero General Public License  *
* along with this program.  If not, see <https://www.gnu.org/licenses/>.    *
*****************************************************************************/

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