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
// NAME:              Fractures.cuh
//
// PURPOSE:           Create fractures in a DFN
//                    ModeSizeDistri, // 0 = power law; 1 = lognormal; 2 = uniform; 3 = monosize
//                    ParaSizeDistri: when 0, ParaSizeDistri.x = alpha, y = minR, z = maxR
//                                    when 1, x = mean, y = sigma, z = minR, w = maxR
//                                    when 2, x = minR, y = minR
//                                    when 3, x = R
//
// FUNCTIONS/OBJECTS: Fractures
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../DataTypeSelector/DataTypeSelector.cuh"
#include "../GlobalDef/GlobalDef.cuh"
#include "../MatrixManipulation/MatrixManipulation.cuh"
#include "../RandomFunction/RandomFunction.cuh"
#include "Fracture.cuh"
#include "TruncateFracture.cuh"

namespace cuDFNsys
{
// Generate some fractures in a DFN
template <typename T>
__global__ void Fractures(cuDFNsys::Fracture<T> *verts,
                          unsigned long seed,
                          int count,
                          cuDFNsys::Vector1<T> model_L,
                          uint ModeSizeDistri,                 // 0 = power law; 1 = lognormal; 2 = uniform; 3 = monosize
                          cuDFNsys::Vector4<T> ParaSizeDistri, // when mode = 1, ;
                          cuDFNsys::Vector1<T> kappa,
                          cuDFNsys::Vector1<T> conductivity_powerlaw_exponent,
                          T Gamma_constant = 5.0e-4);

// Generate some fractures in a DFN
template <typename T>
class FracturesCPU
{
public:
    // constructor
    FracturesCPU(thrust::host_vector<cuDFNsys::Fracture<T>> &verts,
                 unsigned long seed,
                 int count,
                 T model_L,
                 uint ModeSizeDistri,                 // 0 = power law; 1 = lognormal; 2 = uniform; 3 = monosize
                 cuDFNsys::Vector4<T> ParaSizeDistri, // when mode = 1, ;
                 T kappa,
                 T conductivity_powerlaw_exponent,
                 uint Nproc = 10);
};

// benchmark fracture generator: two vertical crossed fractures
template <typename T>
__global__ void FracturesCrossedVertical(cuDFNsys::Fracture<T> *verts,
                                         unsigned long seed,
                                         int count,
                                         T model_L);

// benchmark fracture generator: two inclined fractures, with two beta values
// beta is dip angle here
template <typename T>
__global__ void FracturesBeta50Beta60(cuDFNsys::Fracture<T> *verts,
                                      unsigned long seed,
                                      int count,
                                      T model_L);

// two incomplet fractures
template <typename T>
__global__ void FracturesIncomplete(cuDFNsys::Fracture<T> *verts,
                                    unsigned long seed,
                                    int count,
                                    T model_L);

// Fractures like 2D sticks
template <typename T>
__global__ void Fractures2DLike(cuDFNsys::Fracture<T> *verts,
                                unsigned long seed,
                                int count,
                                T model_L,
                                T alpha = 1.5,
                                T minR = 1,
                                T maxR = 15);

// Four fractures to verify particle tracking
template <typename T>
__global__ void FracturesFour(cuDFNsys::Fracture<T> *verts,
                              unsigned long seed,
                              int count,
                              T model_L);

template <typename T>
__global__ void FracturesChangeDomainSize(cuDFNsys::Fracture<T> *verts,
                                          int count,
                                          T model_L);
}; // namespace cuDFNsys