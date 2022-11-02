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

#include "Geometry/2D/Intersection2DLine2DSeg.cuh"

// ====================================================
// NAME:        Intersection2DLine2DSeg
// DESCRIPTION: Identify Intersections between
//              2D line (infinite) and 2D segment
// AUTHOR:      Tingchang YIN
// DATE:        07/04/2022
// ====================================================
template <typename T>
__device__ __host__ bool cuDFNsys::Intersection2DLine2DSeg(cuDFNsys::Vector2<T> *Line,
                                                           cuDFNsys::Vector2<T> *Seg,
                                                           int *sign_, // 1, pnt; 2, seg; 3, none
                                                           cuDFNsys::Vector2<T> *intersection,
                                                           T _TOL_)
{
    cuDFNsys::Vector2<T> directional_v = cuDFNsys::MakeVector2(Line[0].x - Line[1].x,
                                                               Line[0].y - Line[1].y);

    cuDFNsys::Vector2<T> Normal_To_Line = cuDFNsys::MakeVector2(-1.0f * directional_v.y,
                                                                directional_v.x);

    cuDFNsys::Vector2<T> p1 = Seg[0];
    cuDFNsys::Vector2<T> p2 = Seg[1];

    cuDFNsys::Vector2<T> p = cuDFNsys::MakeVector2((Line[0].x + Line[1].x) * 0.5f,
                                                   (Line[0].y + Line[1].y) * 0.5f);

    cuDFNsys::Vector2<T> j1 = cuDFNsys::MakeVector2(p1.x - p.x, p1.y - p.y);
    cuDFNsys::Vector2<T> j2 = cuDFNsys::MakeVector2(p2.x - p.x, p2.y - p.y);

    T g1 = Normal_To_Line.x * j1.x + Normal_To_Line.y * j1.y;
    T g2 = Normal_To_Line.x * j2.x + Normal_To_Line.y * j2.y;

    //printf("%f %f\n", g1, g2);

    if (abs(g1) < _TOL_ && abs(g2) < _TOL_) // L is containing the Seg
    {
        *sign_ = 2;
        intersection[0] = Seg[0];
        intersection[1] = Seg[1];
        return true;
    }
    else if (abs(g1) < _TOL_ && abs(g2) > _TOL_) // one point lies on the L
    {
        *sign_ = 1;
        intersection[0] = Seg[0];
        return true;
    }
    else if (abs(g1) > _TOL_ && abs(g2) < _TOL_) // another point lies on the L
    {
        *sign_ = 1;
        intersection[0] = Seg[1];
        return true;
    }

    if ((g1 / abs(g1)) + (g2 / abs(g2)) == 0) // oppsite signs
    {
        *sign_ = 1;
        // now intersection between two infinite lines is the intersection we need
        T a1 = Line[0].x;
        T a2 = Line[0].y;
        T b1 = Line[1].x;
        T b2 = Line[1].y;

        T c1 = Seg[0].x;
        T c2 = Seg[0].y;
        T d1 = Seg[1].x;
        T d2 = Seg[1].y;

        T A1 = a2 - b2;
        T B1 = b1 - a1;
        T C1 = (A1 * a1) + (B1 * a2);
        T A2 = c2 - d2;
        T B2 = d1 - c1;
        T C2 = (A2 * c1) + (B2 * c2);

        T x = (C1 * B2 - C2 * B1) / (A1 * B2 - A2 * B1);
        T y = (C2 * A1 - C1 * A2) / (A1 * B2 - B1 * A2);
        intersection[0] = cuDFNsys::MakeVector2(x, y);
        return true;
    }

    *sign_ = -1;
    return false;
}; // Intersection2DLine2DSeg
template __device__ __host__ bool cuDFNsys::Intersection2DLine2DSeg<double>(cuDFNsys::Vector2<double> *Line,
                                                                            cuDFNsys::Vector2<double> *Seg,
                                                                            int *sign_, // 1, pnt; 2, seg; 3, none
                                                                            cuDFNsys::Vector2<double> *intersection,
                                                                            double _TOL_);
template __device__ __host__ bool cuDFNsys::Intersection2DLine2DSeg<float>(cuDFNsys::Vector2<float> *Line,
                                                                           cuDFNsys::Vector2<float> *Seg,
                                                                           int *sign_, // 1, pnt; 2, seg; 3, none
                                                                           cuDFNsys::Vector2<float> *intersection,
                                                                           float _TOL_);