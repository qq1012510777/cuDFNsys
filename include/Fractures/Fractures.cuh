///////////////////////////////////////////////////////////////////
// NAME:              Fractures.cuh
//
// PURPOSE:           Create fractures in a DFN
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
                          float alpha,
                          float minR,
                          float maxR,
                          float kappa = 0,
                          float conductivity_powerlaw_exponent = 0);

}; // namespace cuDFNsys