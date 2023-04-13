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
// NAME:              IdentifyIntersection.cuh
//
// PURPOSE:           Identify intersection
//
// FUNCTIONS/OBJECTS: IdentifyIntersection
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../Exceptions/Exceptions.cuh"
#include "../Geometry/Geometry.cuh"
#include "../GlobalDef/GlobalDef.cuh"
#include "../MatrixManipulation/MatrixManipulation.cuh"
#include "Fracture.cuh"
#include "IdentifyFracPairSphericalDetection.cuh"
#include "IdentifyIntersectionKernel.cuh"
#include "Intersection.cuh"
#include <unistd.h>

namespace cuDFNsys
{
template <typename T>
class IdentifyIntersection
{
public:
    // constructor CPU
    IdentifyIntersection(thrust::host_vector<cuDFNsys::Fracture<T>> verts,
                         const bool &if_trucncated,
                         std::map<pair<size_t, size_t>, pair<cuDFNsys::Vector3<T>, cuDFNsys::Vector3<T>>> &Intersection_map, uint Nproc = 10);
    // constructor GPU
    IdentifyIntersection(const size_t &Fracsize,
                         cuDFNsys::Fracture<T> *Frac_verts_device_ptr,
                         const bool &if_trucncated,
                         std::map<pair<size_t, size_t>, pair<cuDFNsys::Vector3<T>, cuDFNsys::Vector3<T>>> &Intersection_map);
};
}; // namespace cuDFNsys