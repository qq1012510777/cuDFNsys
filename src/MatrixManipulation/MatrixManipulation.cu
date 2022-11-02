/****************************************************************************
* cuDFNsys - simulating flow and transport in 3D fracture networks          *
* Copyright (C) 2022, Tingchang YIN, Sergio GALINDO-TORRES                  *
*                                                                           *
* This program is free software: you can redistribute it and/or modify      *
* it under the terms of the GNU Affero General Public License as            *
* published by the Free Software Foundation, either version 3 of the        *
* License, or (at your option) any later version.                           *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU Affero General Public License for more details.                       *
*                                                                           *
* You should have received a copy of the GNU Affero General Public License  *
* along with this program.  If not, see <https://www.gnu.org/licenses/>.    *
*****************************************************************************/

#include "MatrixManipulation/MatrixManipulation.cuh"
// ====================================================
// NAME:        MakeVector2
// DESCRIPTION: create vector2
// AUTHOR:      Tingchang YIN
// DATE:        02/07/2022
// ====================================================
template <typename T>
__device__ __host__ cuDFNsys::Vector2<T> cuDFNsys::MakeVector2(T a, T b)
{
    cuDFNsys::Vector2<T> result;
    result.x = a, result.y = b;
    return result;
}; // MakeVector2
template __device__ __host__ cuDFNsys::Vector2<double> cuDFNsys::MakeVector2(cuDFNsys::Vector1<double> a, cuDFNsys::Vector1<double> b);
template __device__ __host__ cuDFNsys::Vector2<float> cuDFNsys::MakeVector2(cuDFNsys::Vector1<float> a, cuDFNsys::Vector1<float> b);

// ====================================================
// NAME:        MakeVector3
// DESCRIPTION: create vector3
// AUTHOR:      Tingchang YIN
// DATE:        02/07/2022
// ====================================================
template <typename T>
__device__ __host__ cuDFNsys::Vector3<T> cuDFNsys::MakeVector3(T a,
                                                               T b,
                                                               T c)
{
    cuDFNsys::Vector3<T> result;
    result.x = a, result.y = b, result.z = c;
    return result;
}; // MakeVector3
template __device__ __host__ cuDFNsys::Vector3<double> cuDFNsys::MakeVector3(cuDFNsys::Vector1<double> a, cuDFNsys::Vector1<double> b, cuDFNsys::Vector1<double> c);
template __device__ __host__ cuDFNsys::Vector3<float> cuDFNsys::MakeVector3(cuDFNsys::Vector1<float> a, cuDFNsys::Vector1<float> b, cuDFNsys::Vector1<float> c);

// ====================================================
// NAME:        MakeVector4
// DESCRIPTION: create vector4
// AUTHOR:      Tingchang YIN
// DATE:        02/07/2022
// ====================================================
template <typename T>
__device__ __host__ cuDFNsys::Vector4<T> cuDFNsys::MakeVector4(T a, T b, T c, T d)
{
    cuDFNsys::Vector4<T> result;
    result.x = a, result.y = b, result.z = c, result.w = d;
    return result;
}; // MakeVector4
template __device__ __host__ cuDFNsys::Vector4<double> cuDFNsys::MakeVector4(cuDFNsys::Vector1<double> a, cuDFNsys::Vector1<double> b, cuDFNsys::Vector1<double> c, cuDFNsys::Vector1<double> d);
template __device__ __host__ cuDFNsys::Vector4<float> cuDFNsys::MakeVector4(cuDFNsys::Vector1<float> a, cuDFNsys::Vector1<float> b, cuDFNsys::Vector1<float> c, cuDFNsys::Vector1<float> d);

// ====================================================
// NAME:        CrossProductVector3
// DESCRIPTION: Cross product of two 3D vectors.
// AUTHOR:      Tingchang YIN
// DATE:        04/04/2022
// ====================================================
template <typename T>
__device__ __host__ cuDFNsys::Vector3<T> cuDFNsys::CrossProductVector3(cuDFNsys::Vector3<T> v1,
                                                                       cuDFNsys::Vector3<T> v2)
{
    cuDFNsys::Vector3<T> n;
    n.x = v1.y * v2.z - v1.z * v2.y,
    n.y = v1.z * v2.x - v1.x * v2.z,
    n.z = v1.x * v2.y - v1.y * v2.x;
    return n;
}; // CrossProductVector3
template __device__ __host__ cuDFNsys::Vector3<double> cuDFNsys::CrossProductVector3<double>(cuDFNsys::Vector3<double> v1, cuDFNsys::Vector3<double> v2);
template __device__ __host__ cuDFNsys::Vector3<float> cuDFNsys::CrossProductVector3<float>(cuDFNsys::Vector3<float> v1, cuDFNsys::Vector3<float> v2);

// ====================================================
// NAME:        CroCrossProductVector2
// DESCRIPTION: because the A and B are in XY plane,
//              so the cross product is
//              (0, 0, A.x * B.y - B.x * A.y).
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================
template <typename T>
__device__ __host__ T cuDFNsys::CrossProductVector2(cuDFNsys::Vector2<T> A, cuDFNsys::Vector2<T> B)
{
    return A.x * B.y - B.x * A.y;
}; // CroCrossProductVector2
template __device__ __host__ cuDFNsys::Vector1<double> cuDFNsys::CrossProductVector2<double>(cuDFNsys::Vector2<double> A, cuDFNsys::Vector2<double> B);
template __device__ __host__ cuDFNsys::Vector1<float> cuDFNsys::CrossProductVector2<float>(cuDFNsys::Vector2<float> A, cuDFNsys::Vector2<float> B);

// ====================================================
// NAME:        ProductSquare3Vector3
// DESCRIPTION: Product of a square matrix and a column vector.
// AUTHOR:      Tingchang YIN
// DATE:        07/04/2022
// ====================================================
template <typename T>
__device__ __host__ cuDFNsys::Vector3<T> cuDFNsys::ProductSquare3Vector3(T A[3][3], cuDFNsys::Vector3<T> B)
{
    T B_[3] = {B.x, B.y, B.z};
    T C[3] = {0};

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            C[i] += A[i][j] * B_[j];
        }
    }

    cuDFNsys::Vector3<T> C_ = cuDFNsys::MakeVector3(C[0], C[1], C[2]);

    //C_.x = C[0], C_.y = C[1], C_.z = C[2];

    return C_;
}; // ProductSquare3Vector3
template __device__ __host__ cuDFNsys::Vector3<double> cuDFNsys::ProductSquare3Vector3<double>(double A[3][3],
                                                                                     cuDFNsys::Vector3<double> B);
template __device__ __host__ cuDFNsys::Vector3<float> cuDFNsys::ProductSquare3Vector3<float>(float A[3][3],
                                                                                   cuDFNsys::Vector3<float> B);

// ====================================================
// NAME:        ProjectVToPlaneN
// DESCRIPTION: project a vector V to a plane which has normal of n
// AUTHOR:      Tingchang YIN
// DATE:        12/05/2022
// ====================================================
template <typename T>
__device__ __host__ cuDFNsys::Vector3<T> cuDFNsys::ProjectVToPlaneN(cuDFNsys::Vector3<T> V, cuDFNsys::Vector3<T> n)
{
    T I_minus_nn[3][3] = {1.0f - n.x * n.x, -n.x * n.y, -n.x * n.z,
                          -n.y * n.x, 1.0f - n.y * n.y, -n.y * n.z,
                          -n.z * n.x, -n.z * n.y, 1.0f - n.z * n.z};

    cuDFNsys::Vector3<T> KL = cuDFNsys::ProductSquare3Vector3<T>(I_minus_nn, V);
    T norm_KL = sqrt(KL.x * KL.x + KL.y * KL.y + KL.z * KL.z);

    KL.x /= norm_KL;
    KL.y /= norm_KL;
    KL.z /= norm_KL;

    return KL;
}; // ProjectVToPlaneN
template __device__ __host__ cuDFNsys::Vector3<double> cuDFNsys::ProjectVToPlaneN<double>(cuDFNsys::Vector3<double> V,
                                                                                cuDFNsys::Vector3<double> n);
template __device__ __host__ cuDFNsys::Vector3<float> cuDFNsys::ProjectVToPlaneN<float>(cuDFNsys::Vector3<float> V,
                                                                              cuDFNsys::Vector3<float> n);