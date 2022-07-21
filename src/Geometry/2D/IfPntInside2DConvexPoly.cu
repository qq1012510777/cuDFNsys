#include "Geometry/2D/IfPntInside2DConvexPoly.cuh"
// ====================================================
// NAME:        IfPntInside2DConvexPoly
// DESCRIPTION: if a point is inside a 2D polygon
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================
template <typename T>
__device__ __host__ bool cuDFNsys::IfPntInside2DConvexPoly(cuDFNsys::Vector2<T> pnt,
                                                           cuDFNsys::Vector2<T> *verts,
                                                           int N)
{
    bool inside = true;

    T R[10];

    for (int i = 0; i < N; ++i)
    {
        cuDFNsys::Vector2<T> A = cuDFNsys::MakeVector2(verts[i].x - pnt.x, verts[i].y - pnt.y);
        T norm_A = pow(A.x * A.x + A.y * A.y, 0.5);
        A = cuDFNsys::MakeVector2(A.x / norm_A, A.y / norm_A);

        cuDFNsys::Vector2<T> B = cuDFNsys::MakeVector2(verts[(i + 1) % N].x - pnt.x, verts[(i + 1) % N].y - pnt.y);
        T norm_B = pow(B.x * B.x + B.y * B.y, 0.5);
        B = cuDFNsys::MakeVector2(B.x / norm_B, B.y / norm_B);

        R[i] = cuDFNsys::CrossProductVector2<T>(A, B);

        R[i] = R[i] / abs(R[i]);
    };

    for (int i = 0; i < N; ++i)
    {
        if (R[i] != R[(i + 1) % N])
        {
            inside = false;
            break;
        }
    }

    return inside;
}; // IfPntInside2DConvexPoly
template __device__ __host__ bool cuDFNsys::IfPntInside2DConvexPoly<double>(cuDFNsys::Vector2<double> pnt,
                                                                            cuDFNsys::Vector2<double> *verts,
                                                                            int N);
template __device__ __host__ bool cuDFNsys::IfPntInside2DConvexPoly<float>(cuDFNsys::Vector2<float> pnt,
                                                                           cuDFNsys::Vector2<float> *verts,
                                                                           int N);