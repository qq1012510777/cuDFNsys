///////////////////////////////////////////////////////////////////
// NAME:              GlobalDef.cuh
//
// PURPOSE:           A Warmup funtion
//
// FUNCTIONS/OBJECTS: N/A
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../GlobalDef/GlobalDef.cuh"

// warmup GPU
namespace cuDFNsys
{
__global__ void Warmup();
};