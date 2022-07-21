#include "Geometry/2D/IfPntLiesOnBound2DConvexPoly.cuh"

// ====================================================
// NAME:        IfPntLiesOnBound2DConvexPoly
// DESCRIPTION: if a point lies on 2D bound
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================
template <typename T>
__device__ __host__ bool cuDFNsys::IfPntLiesOnBound2DConvexPoly(cuDFNsys::Vector2<T> pnt,
                                                                cuDFNsys::Vector2<T> *verts,
                                                                int N,
                                                                T _tol_)
{
    for (int i = 0; i < N; ++i)
    {
        cuDFNsys::Vector2<T> Seg[2] = {verts[i], verts[(i + 1) % N]};

        T dist = cuDFNsys::DistancePnt2DSeg<T>(pnt, Seg);

        if (dist < _tol_)
            return true;
    }
    return false;
}; // IfPntLiesOnBound2DConvexPoly
template __device__ __host__ bool cuDFNsys::IfPntLiesOnBound2DConvexPoly<double>(cuDFNsys::Vector2<double> pnt,
                                                                                 cuDFNsys::Vector2<double> *verts,
                                                                                 int N,
                                                                                 double _tol_);
template __device__ __host__ bool cuDFNsys::IfPntLiesOnBound2DConvexPoly<float>(cuDFNsys::Vector2<float> pnt,
                                                                                cuDFNsys::Vector2<float> *verts,
                                                                                int N,
                                                                                float _tol_);

// ====================================================
// NAME:        IfPntLiesOnBound2DConvexPolyReturnEdgeNO
// DESCRIPTION: if a point lies on 2D bound
//              also return the edge NO
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================
template <typename T>
__device__ __host__ bool cuDFNsys::IfPntLiesOnBound2DConvexPolyReturnEdgeNO(cuDFNsys::Vector2<T> pnt,
                                                                            cuDFNsys::Vector2<T> *verts,
                                                                            int N,
                                                                            T _tol_,
                                                                            int *edgeNO) // 0 1 2 3 4 5
{
    for (int i = 0; i < N; ++i)
    {
        cuDFNsys::Vector2<T> Seg[2] = {verts[i], verts[(i + 1) % N]};

        T dist = cuDFNsys::DistancePnt2DSeg<T>(pnt, Seg);

        if (dist < _tol_)
        {
            *edgeNO = i;
            return true;
        }
    }

    return false;
}; // IfPntLiesOnBound2DConvexPolyReturnEdgeNO
template __device__ __host__ bool cuDFNsys::IfPntLiesOnBound2DConvexPolyReturnEdgeNO<double>(cuDFNsys::Vector2<double> pnt,
                                                                                             cuDFNsys::Vector2<double> *verts,
                                                                                             int N,
                                                                                             double _tol_,
                                                                                             int *edgeNO);
template __device__ __host__ bool cuDFNsys::IfPntLiesOnBound2DConvexPolyReturnEdgeNO<float>(cuDFNsys::Vector2<float> pnt,
                                                                                            cuDFNsys::Vector2<float> *verts,
                                                                                            int N,
                                                                                            float _tol_,
                                                                                            int *edgeNO);