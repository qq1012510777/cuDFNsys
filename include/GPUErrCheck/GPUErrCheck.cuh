///////////////////////////////////////////////////////////////////
// NAME:              GPUErrCheck.cuh
//
// PURPOSE:           A GPUErrCheck funtion: check cuda_success
//
// FUNCTIONS/OBJECTS: N/A
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../GlobalDef/GlobalDef.cuh"

#define GPUErrCheck(ans)                      \
    {                                         \
        gpuAssert((ans), __FILE__, __LINE__); \
    }
inline void gpuAssert(cudaError_t code, std::string file, int line, bool abort = true)
{
    if (code != cudaSuccess)
    {
        fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file.c_str(), line);
        if (abort)
            exit(code);
    }
};