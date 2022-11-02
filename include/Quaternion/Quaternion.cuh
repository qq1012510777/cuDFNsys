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
// NAME:              Quaternion.cuh
//
// PURPOSE:           Quaternion function to
//                    rotate a point around an axis
// FUNCTIONS/OBJECTS: N/A
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../DataTypeSelector/DataTypeSelector.cuh"
#include "../GlobalDef/GlobalDef.cuh"

// Quaternion helper class describing rotations
// this allows for a nice description and execution of rotations in 3D space.
namespace cuDFNsys
{
template <typename T>
struct Quaternion
{
protected:
    // 1,i,j,k
    cuDFNsys::Vector4<T> QuaternionNum;

public:
    // describe quaternion
    __device__ __host__ Quaternion DescribeRotation(const cuDFNsys::Vector3<T> v, const cuDFNsys::Vector1<T> angle);

    // rotate
    __device__ __host__ cuDFNsys::Vector3<T> Rotate(const cuDFNsys::Vector3<T> v);

    // get QuaternionNum
    __device__ __host__ cuDFNsys::Vector4<T> GetQuaternionNum();
};
}; // namespace cuDFNsys