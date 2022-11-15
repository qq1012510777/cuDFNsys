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

#include "ParticleTransport/ParticleReflection.cuh"

// ====================================================
// NAME:        ParticleReflection
// DESCRIPTION: apply the specular rule for particles
// AUTHOR:      Tingchang YIN
// DATE:        15/11/2022
// ====================================================

template <typename T>
__device__ __host__ cuDFNsys::Vector2<T> cuDFNsys::ParticleReflection(cuDFNsys::Vector2<T> P,
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
}; // ParticleReflection
template __device__ __host__ cuDFNsys::Vector2<double> cuDFNsys::ParticleReflection<double>(cuDFNsys::Vector2<double> P,
                                                                                            cuDFNsys::Vector2<double> A,
                                                                                            cuDFNsys::Vector2<double> B);
template __device__ __host__ cuDFNsys::Vector2<float> cuDFNsys::ParticleReflection<float>(cuDFNsys::Vector2<float> P,
                                                                                          cuDFNsys::Vector2<float> A,
                                                                                          cuDFNsys::Vector2<float> B);
