///////////////////////////////////////////////////////////////////
// NAME:              Quaternion.cuh
//
// PURPOSE:           Quaternion function to
//                    rotate a point around an axis
// FUNCTIONS/OBJECTS: N/A
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../GlobalDef/GlobalDef.cuh"

// Quaternion helper class describing rotations
// this allows for a nice description and execution of rotations in 3D space.
namespace cuDFNsys
{
struct Quaternion
{
protected:
    // 1,i,j,k
    float4 QuaternionNum;

public:
    __device__ __host__ Quaternion DescribeRotation(const float3 v, const float angle);
    __device__ __host__ float3 Rotate(const float3 v);
    __device__ __host__ float4 GetQuaternionNum();
};
}; // namespace cuDFNsys