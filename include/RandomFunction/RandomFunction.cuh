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
// NAME:              RandomFunction.cuh
//
// PURPOSE:           some functions to generate random numbers
//
// FUNCTIONS/OBJECTS: N/A
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../GlobalDef/GlobalDef.cuh"

namespace cuDFNsys
{
// power law
template <typename T>
__device__ __host__ T RandomPowerlaw(T x0,
                                     T x1,
                                     T alpha_g,
                                     T rand_0_1);

// uniform
template <typename T>
__device__ __host__ T RandomUniform(T a,
                                    T b,
                                    T rand_0_1);

// fisher
template <typename T>
__device__ __host__ T RandomFisher(T rand_0_1,
                                   T fisher_k);

// lognormal: mean is the mean of log(x),
// variance is the variance of log(x)
template <typename T>
__device__ __host__ T RandomLognormal(T mean,
                                      T variance,
                                      T min,
                                      T max,
                                      T rand_0_1);

// return CDF of standard normal distribution (-inf, x)
template <typename T>
__device__ __host__ T StandardNormalCDF(T x);

// return a value that is calculated by inverting the standard normal CDF.
template <typename T>
__device__ __host__ T StandardNormalCDFInv(T p);

// a stable algorithm to calculate polynomial with 8th degrees
template <typename T>
__device__ __host__ T R8PolyValueHorner(int m, T c[], T x);
}; // namespace cuDFNsys