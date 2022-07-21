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
#include "../DataTypeSelector/DataTypeSelector.cuh"
#include "../GlobalDef/GlobalDef.cuh"

// Quaternion helper class describing rotations
// this allows for a nice description and execution of rotations in 3D space.
namespace cuDFNsys
{
template <typename T>
struct Quaternion
{
protected:
    // 1,i,j,k
    cuDFNsys::Vector4<T> QuaternionNum;

public:
    // describe quaternion
    __device__ __host__ Quaternion DescribeRotation(const cuDFNsys::Vector3<T> v, const cuDFNsys::Vector1<T> angle);

    // rotate
    __device__ __host__ cuDFNsys::Vector3<T> Rotate(const cuDFNsys::Vector3<T> v);

    // get QuaternionNum
    __device__ __host__ cuDFNsys::Vector4<T> GetQuaternionNum();
};
}; // namespace cuDFNsys