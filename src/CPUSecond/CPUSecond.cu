#include "CPUSecond/CPUSecond.cuh"
// ====================================================
// NAME:        CPUSecond
// DESCRIPTION: get the current time
//
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================
double cuDFNsys::CPUSecond()
{
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return ((double)tp.tv_sec + (double)tp.tv_usec * 1e-6);
}; // CPUSecond