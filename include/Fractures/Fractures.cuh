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
#include "../GlobalDef/GlobalDef.cuh"
#include "../MatrixManipulation/MatrixManipulation.cuh"
#include "../RandomFunction/RandomFunction.cuh"
#include "Fracture.cuh"
#include "TruncateFracture.cuh"

namespace cuDFNsys
{
// Generate some fractures in a DFN
__global__ void Fractures(cuDFNsys::Fracture *verts,
                          unsigned long seed,
                          int count,
                          float model_L,
                          uint ModeSizeDistri,   // 1 = power law; 2 = lognormal; 3 = uniform; 4 = monosize
                          float4 ParaSizeDistri, // when mode = 1, ;
                          float kappa = 0,
                          float conductivity_powerlaw_exponent = 0);

}; // namespace cuDFNsys