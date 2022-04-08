#include "CpuSecond/CpuSecond.cuh"
// ====================================================
// NAME:        CpuSecond
// DESCRIPTION: get the current time
//
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================
double cuDFNsys::CpuSecond()
{
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return ((double)tp.tv_sec + (double)tp.tv_usec * 1e-6);
}; // CpuSecond