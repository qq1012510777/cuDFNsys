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
        //printf("dist: %.40f, _tol_: %.40f\n", dist, _tol_);
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
        // printf("IfPntLiesOnBound2DConvexPolyReturnEdgeNO: dist: %.40f\n", dist);
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