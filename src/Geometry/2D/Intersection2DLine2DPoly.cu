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

#include "Geometry/2D/Intersection2DLine2DPoly.cuh"

// ====================================================
// NAME:        Intersection2DLine2DPoly
// DESCRIPTION: Identify Intersections between
//              2D Line (infinite) and 2D Polygon
// AUTHOR:      Tingchang YIN
// DATE:        07/04/2022
// ====================================================
template <typename T>
__device__ __host__ bool cuDFNsys::Intersection2DLine2DPoly(cuDFNsys::Vector2<T> *Poly2D,
                                                            int NUM_verts,
                                                            cuDFNsys::Vector2<T> *Line,
                                                            cuDFNsys::Vector2<T> *intersection_k,
                                                            T _TOL_)
{
    int tmpcc = 0;
    cuDFNsys::Vector2<T> tmp_pnts[10];
    for (int i = 0; i < NUM_verts; ++i)
    {
        cuDFNsys::Vector2<T> Seg[2], intersection[2];
        Seg[0] = Poly2D[i];
        Seg[1] = Poly2D[(i + 1) % NUM_verts];

        int sign__ = 0;

        bool ik = cuDFNsys::Intersection2DLine2DSeg<T>(Line,
                                                       Seg,
                                                       &sign__,
                                                       intersection,
                                                       _TOL_);

        if (ik == true)
        {

            if (sign__ == 1)
            {

                bool if_duplicated = false;
                for (int j = 0; j < tmpcc; ++j)
                {
                    cuDFNsys::Vector2<T> dd = cuDFNsys::MakeVector2(intersection[0].x - tmp_pnts[j].x, intersection[0].y - tmp_pnts[j].y);
                    T distance = pow(dd.x * dd.x + dd.y * dd.y, 0.5);

                    if (distance < _TOL_)
                    {
                        if_duplicated = true;
                        break;
                    }
                }

                if (if_duplicated == false)
                {
                    tmp_pnts[tmpcc] = intersection[0];
                    tmpcc++;
                }
            }
        }
    }

    if (tmpcc == 2)
    {
        intersection_k[0] = tmp_pnts[0];
        intersection_k[1] = tmp_pnts[1];
        return true;
    }

    if (tmpcc == 1)
    {
        return false;
    }

    if (tmpcc == 0)
    {
        return false;
    }

    if (tmpcc > 2)
    {
        //printf("Error! An infinite line (2D) cannot intersect a 2D polygon with more than two points!\n");
        //exit(0);

        uint2 pair__;

        T dist = 0;
        for (int k = 0; k < tmpcc - 1; ++k)
        {
            for (int h = k + 1; h < tmpcc; ++h)
            {
                cuDFNsys::Vector2<T> kh = cuDFNsys::MakeVector2(tmp_pnts[k].x - tmp_pnts[h].x, tmp_pnts[k].y - tmp_pnts[h].y);
                T distrr = sqrt(kh.x * kh.x + kh.y * kh.y);

                if (distrr > dist)
                {
                    dist = distrr;
                    pair__.x = k;
                    pair__.y = h;
                }
            }
        };

        intersection_k[0] = tmp_pnts[pair__.x];
        intersection_k[1] = tmp_pnts[pair__.y];
        return true;
    }

    return false;
}; // Intersection2DLine2DPoly
template __device__ __host__ bool cuDFNsys::Intersection2DLine2DPoly<double>(cuDFNsys::Vector2<double> *Poly2D,
                                                                    int NUM_verts,
                                                                    cuDFNsys::Vector2<double> *Line,
                                                                    cuDFNsys::Vector2<double> *intersection_k,
                                                                    double _TOL_);
template __device__ __host__ bool cuDFNsys::Intersection2DLine2DPoly<float>(cuDFNsys::Vector2<float> *Poly2D,
                                                                   int NUM_verts,
                                                                   cuDFNsys::Vector2<float> *Line,
                                                                   cuDFNsys::Vector2<float> *intersection_k,
                                                                   float _TOL_);
