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
template <typename T>
__device__ __host__ cuDFNsys::Quaternion<T> cuDFNsys::Quaternion<T>::DescribeRotation(const cuDFNsys::Vector3<T> v, const cuDFNsys::Vector1<T> angle)
{
    cuDFNsys::Vector1<T> sina_2 = sin(angle * 0.5000);
    cuDFNsys::Vector1<T> cosa_2 = cos(angle * 0.5000);
    cuDFNsys::Quaternion<T> result;
    result.QuaternionNum.x = cosa_2,
    result.QuaternionNum.y = sina_2 * v.x,
    result.QuaternionNum.z = sina_2 * v.y,
    result.QuaternionNum.w = sina_2 * v.z;

    cuDFNsys::Vector1<T> norm_r = pow(result.QuaternionNum.x * result.QuaternionNum.x + result.QuaternionNum.y * result.QuaternionNum.y +
                                          result.QuaternionNum.z * result.QuaternionNum.z + result.QuaternionNum.w * result.QuaternionNum.w,
                                      0.5);
    result.QuaternionNum.x /= norm_r;
    result.QuaternionNum.y /= norm_r;
    result.QuaternionNum.z /= norm_r;
    result.QuaternionNum.w /= norm_r;
    return result;
}; // Quaternion::DescribeRotation
template __device__ __host__ cuDFNsys::Quaternion<double> cuDFNsys::Quaternion<double>::DescribeRotation(const cuDFNsys::Vector3<double> v, const cuDFNsys::Vector1<double> angle);
template __device__ __host__ cuDFNsys::Quaternion<float> cuDFNsys::Quaternion<float>::DescribeRotation(const cuDFNsys::Vector3<float> v, const cuDFNsys::Vector1<float> angle);

// ====================================================
// NAME:        Rotate
// DESCRIPTION: Rotate a point v around the axis
// AUTHOR:      Tingchang YIN
// DATE:        02/04/2022
// ====================================================
template <typename T>
__device__ __host__ cuDFNsys::Vector3<T> cuDFNsys::Quaternion<T>::Rotate(const cuDFNsys::Vector3<T> v)
{
    cuDFNsys::Vector1<T> t2 = QuaternionNum.x * QuaternionNum.y,
                         t3 = QuaternionNum.x * QuaternionNum.z,
                         t4 = QuaternionNum.x * QuaternionNum.w,
                         t5 = -QuaternionNum.y * QuaternionNum.y,
                         t6 = QuaternionNum.y * QuaternionNum.z,
                         t7 = QuaternionNum.y * QuaternionNum.w,
                         t8 = -QuaternionNum.z * QuaternionNum.z,
                         t9 = QuaternionNum.z * QuaternionNum.w,
                         t10 = -QuaternionNum.w * QuaternionNum.w;

    cuDFNsys::Vector3<T> DF;

    DF.x = 2.0 * ((t8 + t10) * v.x + (t6 - t4) * v.y + (t3 + t7) * v.z) + v.x,
    DF.y = 2.0 * ((t4 + t6) * v.x + (t5 + t10) * v.y + (t9 - t2) * v.z) + v.y,
    DF.z = 2.0 * ((t7 - t3) * v.x + (t2 + t9) * v.y + (t5 + t8) * v.z) + v.z;

    return DF;
}; // Quaternion::Rotate
template __device__ __host__ cuDFNsys::Vector3<double> cuDFNsys::Quaternion<double>::Rotate(const cuDFNsys::Vector3<double> v);
template __device__ __host__ cuDFNsys::Vector3<float> cuDFNsys::Quaternion<float>::Rotate(const cuDFNsys::Vector3<float> v);

// ====================================================
// NAME:        GetQuaternionNum
// DESCRIPTION: get the quaternion number of this struct
// AUTHOR:      Tingchang YIN
// DATE:        02/04/2022
// ====================================================
template <typename T>
__device__ __host__ cuDFNsys::Vector4<T> cuDFNsys::Quaternion<T>::GetQuaternionNum()
{
    cuDFNsys::Vector4<T> f = QuaternionNum;
    return f;
}; // Quaternion::GetQuaternionNum
template __device__ __host__ cuDFNsys::Vector4<double> cuDFNsys::Quaternion<double>::GetQuaternionNum();
template __device__ __host__ cuDFNsys::Vector4<float> cuDFNsys::Quaternion<float>::GetQuaternionNum();