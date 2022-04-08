#include "Quaternion/Quaternion.cuh"

// ====================================================
// NAME:        DescribeRotation
// DESCRIPTION: Generation of Quaternion NUM. 
//              The v should be normalized.
//              v is the rotation axis
//              Angle should be radian
// AUTHOR:      Tingchang YIN
// DATE:        02/04/2022
// ====================================================
__device__ __host__ cuDFNsys::Quaternion cuDFNsys::Quaternion::DescribeRotation(const float3 v, const float angle)
{
    float sina_2 = sin(angle * 0.5000);
    float cosa_2 = cos(angle * 0.5000);
    cuDFNsys::Quaternion result;
    result.QuaternionNum = make_float4(cosa_2, sina_2 * v.x, sina_2 * v.y, sina_2 * v.z);

    float norm_r = pow(result.QuaternionNum.x * result.QuaternionNum.x + result.QuaternionNum.y * result.QuaternionNum.y +
                           result.QuaternionNum.z * result.QuaternionNum.z + result.QuaternionNum.w * result.QuaternionNum.w,
                       0.5);
    result.QuaternionNum.x /= norm_r;
    result.QuaternionNum.y /= norm_r;
    result.QuaternionNum.z /= norm_r;
    result.QuaternionNum.w /= norm_r;
    return result;
}; // Quaternion::DescribeRotation

// ====================================================
// NAME:        Rotate
// DESCRIPTION: Rotate a point v around the axis
// AUTHOR:      Tingchang YIN
// DATE:        02/04/2022
// ====================================================
__device__ __host__ float3 cuDFNsys::Quaternion::Rotate(const float3 v)
{
    float t2 = QuaternionNum.x * QuaternionNum.y;
    float t3 = QuaternionNum.x * QuaternionNum.z;
    float t4 = QuaternionNum.x * QuaternionNum.w;
    float t5 = -QuaternionNum.y * QuaternionNum.y;
    float t6 = QuaternionNum.y * QuaternionNum.z;
    float t7 = QuaternionNum.y * QuaternionNum.w;
    float t8 = -QuaternionNum.z * QuaternionNum.z;
    float t9 = QuaternionNum.z * QuaternionNum.w;
    float t10 = -QuaternionNum.w * QuaternionNum.w;

    return make_float3(
        2.0 * ((t8 + t10) * v.x + (t6 - t4) * v.y + (t3 + t7) * v.z) + v.x,
        2.0 * ((t4 + t6) * v.x + (t5 + t10) * v.y + (t9 - t2) * v.z) + v.y,
        2.0 * ((t7 - t3) * v.x + (t2 + t9) * v.y + (t5 + t8) * v.z) + v.z);
}; // Quaternion::Rotate

// ====================================================
// NAME:        GetQuaternionNum
// DESCRIPTION: get the quaternion number of this struct
// AUTHOR:      Tingchang YIN
// DATE:        02/04/2022
// ====================================================
__device__ __host__ float4 cuDFNsys::Quaternion::GetQuaternionNum()
{
    float4 f = QuaternionNum;
    return f;
}; // Quaternion::GetQuaternionNum