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

#include "ParticleTransport/RandomWalkerGeometryFunctions.cuh"

// ====================================================
// NAME:        IfRandomWalkLieOnBorderTriangle
// DESCRIPTION: IfRandomWalkLieOnBorderTriangle
// AUTHOR:      Tingchang YIN
// DATE:        29/06/2023
// ====================================================
template <typename T>
__host__ __device__ uint2 cuDFNsys::IfRandomWalkLieOnBorderTriangle<T>(cuDFNsys::Vector2<T> p,
                                                                       cuDFNsys::Vector2<T> Triangle[3],
                                                                       T Tol)
{
    uint2 result;
    result.x = 0;

    for (uint i = 0; i < 3; ++i)
    {
        T norm_ab = sqrt((Triangle[i].x - Triangle[(i + 1) % 3].x) *
                             (Triangle[i].x - Triangle[(i + 1) % 3].x) +
                         (Triangle[i].y - Triangle[(i + 1) % 3].y) *
                             (Triangle[i].y - Triangle[(i + 1) % 3].y));
        T norm_ap = sqrt((Triangle[i].x - p.x) *
                             (Triangle[i].x - p.x) +
                         (Triangle[i].y - p.y) *
                             (Triangle[i].y - p.y));
        T norm_bp = sqrt((Triangle[(i + 1) % 3].x - p.x) *
                             (Triangle[(i + 1) % 3].x - p.x) +
                         (Triangle[(i + 1) % 3].y - p.y) *
                             (Triangle[(i + 1) % 3].y - p.y));
        // printf("--s: %.40f\n", abs(norm_ap + norm_bp - norm_ab));
        if (abs(norm_ap + norm_bp - norm_ab) < Tol)
        {
            result.x = 1;
            result.y = i;
            break;
        }
    }

    return result;
}; // IfRandomWalkLieOnBorderTriangle
template __host__ __device__ uint2 cuDFNsys::IfRandomWalkLieOnBorderTriangle<double>(cuDFNsys::Vector2<double> p,
                                                                                     cuDFNsys::Vector2<double> Triangle[3],
                                                                                     double Tol);
template __host__ __device__ uint2 cuDFNsys::IfRandomWalkLieOnBorderTriangle<float>(cuDFNsys::Vector2<float> p,
                                                                                    cuDFNsys::Vector2<float> Triangle[3],
                                                                                    float Tol);

// ====================================================
// NAME:        Sgn
// DESCRIPTION: sign
// AUTHOR:      Tingchang YIN
// DATE:        29/06/2023
// ====================================================
template <typename T>
__host__ __device__ int cuDFNsys::Sgn<T>(const T &x)
{
    return x >= 0 ? x ? 1 : 0 : -1;
}; // Sgn
template __host__ __device__ int cuDFNsys::Sgn<double>(const double &x);
template __host__ __device__ int cuDFNsys::Sgn<float>(const float &x);

// ====================================================
// NAME:        WhichEdgeDoesTrajectoryIntersect
// DESCRIPTION: WhichEdgeDoesTrajectoryIntersect
// AUTHOR:      Tingchang YIN
// DATE:        29/06/2023
// ====================================================
template <typename T>
__host__ __device__ bool cuDFNsys::WhichEdgeDoesTrajectoryIntersect<T>(cuDFNsys::Vector2<T> pp[2],
                                                                       cuDFNsys::Vector2<T> Triangle[3],
                                                                       uint Result[4])
{
    Result[0] = 0;

    for (uint i = 0; i < 3; ++i)
    {
        cuDFNsys::Vector2<T> a_c, b_c, d_c;
        a_c.x = pp[0].x - Triangle[i].x;
        a_c.y = pp[0].y - Triangle[i].y;

        b_c.x = pp[1].x - Triangle[i].x;
        b_c.y = pp[1].y - Triangle[i].y;

        d_c.x = Triangle[(i + 1) % 3].x - Triangle[i].x;
        d_c.y = Triangle[(i + 1) % 3].y - Triangle[i].y;

        cuDFNsys::Vector2<T> b_a, c_a, d_a;
        b_a.x = pp[1].x - pp[0].x;
        b_a.y = pp[1].y - pp[0].y;

        c_a.x = Triangle[i].x - pp[0].x;
        c_a.y = Triangle[i].y - pp[0].y;

        d_a.x = Triangle[(i + 1) % 3].x - pp[0].x;
        d_a.y = Triangle[(i + 1) % 3].y - pp[0].y;

        //if (cuDFNsys::Sgn<T>(a.cross(b, c)) != cuDFNsys::Sgn<T>(a.cross(b, d)) &&
        //    cuDFNsys::Sgn<T>(c.cross(d, a)) != cuDFNsys::Sgn<T>(c.cross(d, b)))
        //
        //{
        //};
        //printf("KKK:\n\t%.40f\n\n\t%.40f\n\n\t%.40f\n\n\t%.40f\n", d_c.x, a_c.y, d_c.y, a_c.x);
        //
        //printf("i: %d, [a, b]:\n\t[%.40f, %.40f]\n\t[%.40f, %.40f]\n", i, pp[0].x, pp[0].y,
        //       pp[1].x, pp[1].y);
        //printf("i: %d, [c, d]:\n\t[%.40f, %.40f]\n\t[%.40f, %.40f]\n", i, Triangle[i].x, Triangle[i].y,
        //       Triangle[(i + 1) % 3].x, Triangle[(i + 1) % 3].y);
        //printf("%.40f and %.40f\n", (cuDFNsys::CrossProductVector2<T>(b_a, c_a)), (cuDFNsys::CrossProductVector2<T>(b_a, d_a)));
        //printf("%.40f and %.40f\n---\n\n", (cuDFNsys::CrossProductVector2<T>(d_c, a_c)), (cuDFNsys::CrossProductVector2<T>(d_c, b_c)));

        if (cuDFNsys::Sgn<T>(cuDFNsys::CrossProductVector2<T>(b_a, c_a)) != cuDFNsys::Sgn<T>(cuDFNsys::CrossProductVector2<T>(b_a, d_a)) &&
            cuDFNsys::Sgn<T>(cuDFNsys::CrossProductVector2<T>(d_c, a_c)) != cuDFNsys::Sgn<T>(cuDFNsys::CrossProductVector2<T>(d_c, b_c)))

        {
            Result[0]++;
            Result[Result[0]] = i;
        };
    }
    // if (cuDFNsys::CrossProductVector2())
    if (Result[0] > 0)
        return true;
    return false;
}; // WhichEdgeDoesTrajectoryIntersect
template __host__ __device__ bool cuDFNsys::WhichEdgeDoesTrajectoryIntersect<double>(cuDFNsys::Vector2<double> pp[2],
                                                                                     cuDFNsys::Vector2<double> Triangle[3],
                                                                                     uint Result[4]);
template __host__ __device__ bool cuDFNsys::WhichEdgeDoesTrajectoryIntersect<float>(cuDFNsys::Vector2<float> pp[2],
                                                                                    cuDFNsys::Vector2<float> Triangle[3],
                                                                                    uint Result[4]);

// ====================================================
// NAME:        IntersectionLineLine2D
// DESCRIPTION: IntersectionLineLine2D
// AUTHOR:      Tingchang YIN
// DATE:        29/06/2023
// ====================================================
template <typename T>
__host__ __device__ cuDFNsys::Vector3<T> cuDFNsys::IntersectionLineLine2D<T>(cuDFNsys::Vector2<T> pp[2],
                                                                             cuDFNsys::Vector2<T> cc[2])
{
    T a1 = pp[1].y - pp[0].y;
    T b1 = pp[0].x - pp[1].x;
    T c1 = a1 * (pp[0].x) + b1 * (pp[0].y);

    // Line CD represented as a2x + b2y = c2
    T a2 = cc[1].y - cc[0].y;
    T b2 = cc[0].x - cc[1].x;
    T c2 = a2 * (cc[0].x) + b2 * (cc[0].y);

    T determinant = a1 * b2 - a2 * b1;

    if (determinant == 0)
    {
        // The lines are parallel. This is simplified
        // by returning a pair of FLT_MAX
        return cuDFNsys::MakeVector3((T)0., (T)0., (T)0.);
    }
    else
    {
        T x = (b2 * c1 - b1 * c2) / determinant;
        T y = (a1 * c2 - a2 * c1) / determinant;
        return cuDFNsys::MakeVector3((T)1, x, y);
    }
}; // IntersectionLineLine2D
template __host__ __device__ cuDFNsys::Vector3<double> cuDFNsys::IntersectionLineLine2D<double>(cuDFNsys::Vector2<double> pp[2],
                                                                                                cuDFNsys::Vector2<double> cc[2]);
template __host__ __device__ cuDFNsys::Vector3<float> cuDFNsys::IntersectionLineLine2D<float>(cuDFNsys::Vector2<float> pp[2],
                                                                                              cuDFNsys::Vector2<float> cc[2]);
