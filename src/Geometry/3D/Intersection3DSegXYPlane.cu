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

#include "Geometry/3D/Intersection3DSegXYPlane.cuh"

// ====================================================
// NAME:        Intersection3DSegXYPlane
// DESCRIPTION: Identify intersection between
//              a 3D segment and the XY plane
// AUTHOR:      Tingchang YIN
// DATE:        07/04/2022
// ====================================================
template <typename T>
__device__ __host__ bool cuDFNsys::Intersection3DSegXYPlane(cuDFNsys::Vector3<T> *Seg,
                                                            cuDFNsys::Vector3<T> *Intersec_PNT,
                                                            int *sign_, // sign: 1: pnt; 2: seg; -1: none;
                                                            T _TOL_)
{
    if (abs(Seg[0].z) < _TOL_ && abs(Seg[1].z) < _TOL_)
    {
        *sign_ = 2;
        Intersec_PNT[0] = Seg[0];
        Intersec_PNT[1] = Seg[1];
        return true;
    }

    if (abs(Seg[0].z) < _TOL_)
    {
        Intersec_PNT[0] = Seg[0];
        *sign_ = 1;
        return true;
    }
    else if (abs(Seg[1].z) < _TOL_)
    {
        Intersec_PNT[0] = Seg[1];
        *sign_ = 1;
        return true;
    }

    T max_z = 0;
    T min_z = 0;

    if (Seg[0].z <= Seg[1].z)
    {
        max_z = Seg[1].z;
        min_z = Seg[0].z;
    }
    else
    {
        min_z = Seg[1].z;
        max_z = Seg[0].z;
    }

    // common cases: cross the xy plane

    if (max_z > 0 && min_z < 0)
    {
        // cross the xy plane
        cuDFNsys::Vector3<T> center_seg = cuDFNsys::MakeVector3((T)(0.5 * (Seg[0].x + Seg[1].x)),
                                                                (T)(0.5 * (Seg[0].y + Seg[1].y)),
                                                                (T)(0.5 * (Seg[0].z + Seg[1].z)));

        cuDFNsys::Vector3<T> vector_seg = cuDFNsys::MakeVector3((T)(Seg[0].x - Seg[1].x),
                                                                (T)(Seg[0].y - Seg[1].y),
                                                                (T)(Seg[0].z - Seg[1].z));

        T vpt = vector_seg.z * 1;

        T t = (0 - center_seg.z) * 1 / vpt;

        cuDFNsys::Vector3<T> Pnt;
        Pnt.x = center_seg.x + vector_seg.x * t;
        Pnt.y = center_seg.y + vector_seg.y * t;
        Pnt.z = center_seg.z + vector_seg.z * t;

        Intersec_PNT[0] = Pnt;

        *sign_ = 1;
        return true;
    }

    *sign_ = -1;
    return false;
}; // Intersection3DSegXYPlane
template __device__ __host__ bool cuDFNsys::Intersection3DSegXYPlane<double>(cuDFNsys::Vector3<double> *Seg,
                                                                             cuDFNsys::Vector3<double> *Intersec_PNT,
                                                                             int *sign_, // sign: 1: pnt; 2: seg; -1: none;
                                                                             double _TOL_);
template __device__ __host__ bool cuDFNsys::Intersection3DSegXYPlane<float>(cuDFNsys::Vector3<float> *Seg,
                                                                            cuDFNsys::Vector3<float> *Intersec_PNT,
                                                                            int *sign_, // sign: 1: pnt; 2: seg; -1: none;
                                                                            float _TOL_);
