#include "MatrixManipulation/MatrixManipulation.cuh"

// ====================================================
// NAME:        CrossProductFloat3
// DESCRIPTION: Cross product of two 3D vectors.
// AUTHOR:      Tingchang YIN
// DATE:        04/04/2022
// ====================================================
__device__ __host__ float3 cuDFNsys::CrossProductFloat3(float3 v1, float3 v2)
{
    float3 n = make_float3(v1.y * v2.z - v1.z * v2.y,
                           v1.z * v2.x - v1.x * v2.z,
                           v1.x * v2.y - v1.y * v2.x);
    return n;
}; // CrossProductFloat3

// ====================================================
// NAME:        CrossProductFloat2
// DESCRIPTION: because the A and B are in XY plane, 
//              so the cross product is 
//              (0, 0, A.x * B.y - B.x * A.y).
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================
__device__ __host__  float cuDFNsys::CrossProductFloat2(float2 A, float2 B)
{
    return A.x * B.y - B.x * A.y;
}

// ====================================================
// NAME:        ProductSquare3Float3
// DESCRIPTION: Product of a square matrix and a column vector.
// AUTHOR:      Tingchang YIN
// DATE:        07/04/2022
// ====================================================
__device__ __host__ float3 cuDFNsys::ProductSquare3Float3(float A[3][3], float3 B)
{
    float B_[3] = {B.x, B.y, B.z};
    float C[3] = {0};

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            C[i] += A[i][j] * B_[j];
        }
    }

    float3 C_ = make_float3(C[0], C[1], C[2]);

    return C_;
}; // ProductSquare3Float3

