#include "Warmup/Warmup.cuh"

// ====================================================
// NAME:        Warmup
// DESCRIPTION: Warmup the GPU.
// AUTHOR:      Tingchang YIN
// DATE:        02/04/2022
// ====================================================
__global__ void cuDFNsys::Warmup()
{
    unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
    float ia, ib;
    ia = ib = 0.0f;
    ib += ia + tid;
}; // Warmup