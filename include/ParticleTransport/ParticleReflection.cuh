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
// NAME:              ParticleReflection.cuh
//
// PURPOSE:           Reflect a particle if the particle crosses a non-flux
//                    edge
//
// FUNCTIONS/OBJECTS: ParticleReflection
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../DataTypeSelector/DataTypeSelector.cuh"
#include "../GlobalDef/GlobalDef.cuh"

namespace cuDFNsys
{
template <typename T>
__device__ __host__ cuDFNsys::Vector2<T> ParticleReflection(cuDFNsys::Vector2<T> P,
                                                            cuDFNsys::Vector2<T> A,
                                                            cuDFNsys::Vector2<T> B)
{
    // Performing translation and shifting origin at A
    cuDFNsys::Vector2<T> Pt = cuDFNsys::MakeVector2(P.x - A.x, P.y - A.y);
    cuDFNsys::Vector2<T> Bt = cuDFNsys::MakeVector2(B.x - A.x, B.y - A.y);

    // Performing rotation in clockwise direction
    // BtAt becomes the X-Axis in the new coordinate system
    cuDFNsys::Vector2<T> Pr = cuDFNsys::MakeVector2((Pt.x * Bt.x + Pt.y * Bt.y) / (Bt.x * Bt.x + Bt.y * Bt.y),
                                                    (T)-1.0 * (Pt.y * Bt.x - Pt.x * Bt.y) / (Bt.x * Bt.x + Bt.y * Bt.y));

    //printf("Pr: %f %f\n", Pr.x, Pr.y);
    // Reflection of Pr about the new X-Axis
    // Followed by restoring from rotation
    // Followed by restoring from translation

    cuDFNsys::Vector2<T> Ps;
    Ps.x = Pr.x * Bt.x - Pr.y * Bt.y + A.x;
    Ps.y = Pr.x * Bt.y + Pr.y * Bt.x + A.y;

    return Ps;
};
template __device__ __host__ cuDFNsys::Vector2<double> ParticleReflection<double>(cuDFNsys::Vector2<double> P,
                                                                                  cuDFNsys::Vector2<double> A,
                                                                                  cuDFNsys::Vector2<double> B);
template __device__ __host__ cuDFNsys::Vector2<float> ParticleReflection<float>(cuDFNsys::Vector2<float> P,
                                                                                cuDFNsys::Vector2<float> A,
                                                                                cuDFNsys::Vector2<float> B);
}; // namespace cuDFNsys