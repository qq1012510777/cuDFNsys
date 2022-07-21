///////////////////////////////////////////////////////////////////
// NAME:              MatrixManipulation.cuh
//
// PURPOSE:           Matrix manipulation
//
// FUNCTIONS/OBJECTS: N/A
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../DataTypeSelector/DataTypeSelector.cuh"
#include "../GlobalDef/GlobalDef.cuh"

namespace cuDFNsys
{
// create cuDFNsys::Vector2<T>
template <typename T>
__device__ __host__ cuDFNsys::Vector2<T> MakeVector2(T a, T b);

// create cuDFNsys::Vector3<T>
template <typename T>
__device__ __host__ cuDFNsys::Vector3<T> MakeVector3(T a,
                                                     T b,
                                                     T c);

// create cuDFNsys::Vector4<T>
template <typename T>
__device__ __host__ cuDFNsys::Vector4<T> MakeVector4(T a, T b, T c, T d);

// cross product of two float3 variables
template <typename T>
__device__ __host__ cuDFNsys::Vector3<T> CrossProductVector3(cuDFNsys::Vector3<T> v1,
                                                             cuDFNsys::Vector3<T> v2);

// cross product of two float2 variables
template <typename T>
__device__ __host__ T CrossProductVector2(cuDFNsys::Vector2<T> A, cuDFNsys::Vector2<T> B);

// product of a square matrix and a column vector
template <typename T>
__device__ __host__ cuDFNsys::Vector3<T> ProductSquare3Vector3(T A[3][3], cuDFNsys::Vector3<T> B);

// project a vector V to a plane which has normal of n
template <typename T>
__device__ __host__ cuDFNsys::Vector3<T> ProjectVToPlaneN(cuDFNsys::Vector3<T> V, cuDFNsys::Vector3<T> n);

}; // namespace cuDFNsys