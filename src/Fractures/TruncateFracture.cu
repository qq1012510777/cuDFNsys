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

#include "Fractures/TruncateFracture.cuh"

// ====================================================
// NAME:        TruncateFracture
// DESCRIPTION: Truncate a fracture in a DFN
// AUTHOR:      Tingchang YIN
// DATE:        06/04/2022
// ====================================================
template <typename T>
__device__ __host__ bool cuDFNsys::TruncateFracture(cuDFNsys::Fracture<T> *verts,
                                                    cuDFNsys::Vector1<T> L,
                                                    int plane,
                                                    int dir)
{
    bool if_touch = false;

    cuDFNsys::Vector1<T> bound = L * 0.5 * dir;
    //printf("bound %f\n", bound);

    cuDFNsys::Vector3<T> TT[8];

    int tmpcc = 0;

    for (int i = 0; i < verts->NumVertsTruncated; ++i)
    {
        cuDFNsys::Vector3<T> source = verts->Verts3DTruncated[i];
        cuDFNsys::Vector3<T> target = verts->Verts3DTruncated[(i + 1) % verts->NumVertsTruncated];

        //printf("The %d edge is being clipped:\n", i + 1);
        //printf("\tsource: %f, %f, %f\n", source.x, source.y, source.z);
        //printf("\ttarget: %f, %f, %f\n", target.x, target.y, target.z);

        cuDFNsys::Vector1<T> x1 = 0;
        cuDFNsys::Vector1<T> x2 = 0;
        if (plane == 0)
        {
            x1 = source.x;
            x2 = target.x;
        }
        else if (plane == 1)
        {
            x1 = source.y;
            x2 = target.y;
        }
        else if (plane == 2)
        {
            x1 = source.z;
            x2 = target.z;
        }

        cuDFNsys::Vector3<T> KK[2];
        int tmpll = -1;

        bool if_both_inside = false;
        bool if_both_outside = false;

        if (dir == -1) // sign = -1
        {
            if (x1 >= bound && x2 >= bound)
            {
                if_both_inside = true;
                if (x1 == bound || x2 == bound)
                    if_touch = true;
            }

            if (x1 < bound && x2 < bound)
                if_both_outside = true;
        }
        else if (dir == 1)
        {
            if (x1 <= bound && x2 <= bound)
            {
                if_both_inside = true;
                if (x1 == bound || x2 == bound)
                    if_touch = true;
            }

            if (x1 > bound && x2 > bound)
                if_both_outside = true;
        }

        //printf("both inside and both outside? %d %d\n", if_both_inside, if_both_outside);

        if (if_both_inside == true)
        {
            KK[0] = source;
            KK[1] = target;
            tmpll = 2;
        }
        else if (if_both_outside == true)
        {
            tmpll = 0;
        }
        else
        {
            cuDFNsys::Vector3<T> discardPNT, keepPNT;

            if (dir == -1) // sign = -1
            {
                if (x1 < bound && x2 >= bound)
                {
                    discardPNT = source;
                    keepPNT = target;
                }
                else if (x1 >= bound && x2 < bound)
                {
                    keepPNT = source;
                    discardPNT = target;
                }
            }
            else if (dir == 1)
            {
                if (x1 > bound && x2 <= bound)
                {
                    discardPNT = source;
                    keepPNT = target;
                }
                else if (x1 <= bound && x2 > bound)
                {
                    keepPNT = source;
                    discardPNT = target;
                }
            }

            cuDFNsys::Vector3<T> v;
            v.x = (discardPNT.x - keepPNT.x),
            v.y = (discardPNT.y - keepPNT.y),
            v.z = (discardPNT.z - keepPNT.z);

            cuDFNsys::Vector1<T> t = 0;
            cuDFNsys::Vector3<T> Intersection_;
            Intersection_.x = 0, Intersection_.y = 0, Intersection_.z = 0;

            if (plane == 0)
            {
                t = (bound - keepPNT.x) / v.x;
                Intersection_.x = bound;
                Intersection_.y = keepPNT.y + t * v.y;
                Intersection_.z = keepPNT.z + t * v.z;
            }
            else if (plane == 1)
            {
                t = (bound - keepPNT.y) / v.y;
                Intersection_.x = keepPNT.x + t * v.x;
                Intersection_.y = bound;
                Intersection_.z = keepPNT.z + t * v.z;
            }
            else if (plane == 2)
            {
                t = (bound - keepPNT.z) / v.z;
                Intersection_.x = keepPNT.x + t * v.x;
                Intersection_.y = keepPNT.y + t * v.y;
                Intersection_.z = bound;
            }

            if (dir == -1)
            {
                if (x1 < bound)
                {
                    KK[0] = Intersection_;
                    KK[1] = keepPNT;
                }
                else if (x2 < bound)
                {
                    KK[1] = Intersection_;
                    KK[0] = keepPNT;
                }
            }
            else if (dir == 1)
            {
                if (x1 > bound)
                {
                    KK[0] = Intersection_;
                    KK[1] = keepPNT;
                }
                else if (x2 > bound)
                {
                    KK[1] = Intersection_;
                    KK[0] = keepPNT;
                }
            }

            tmpll = 2;
            if_touch = true;
        }

        if (tmpll != 0)
        {
            //printf("inserting pnts\n");
            for (int j = 0; j < tmpll; ++j)
            {
                cuDFNsys::Vector3<T> PNT_this = KK[j];

                bool if_duplicated = false;

                for (int k = 0; k < tmpcc; ++k)
                {
                    cuDFNsys::Vector3<T> PNT_e = TT[k];

                    cuDFNsys::Vector3<T> dd;
                    dd.x = PNT_this.x - PNT_e.x, dd.y = PNT_this.y - PNT_e.y, dd.z = PNT_this.z - PNT_e.z;

                    cuDFNsys::Vector1<T> distance = pow(dd.x * dd.x + dd.y * dd.y + dd.z * dd.z, 0.5);
                    if (distance == 0)
                    {
                        if_duplicated = true;
                        break;
                    }
                }

                if (if_duplicated == false)
                {
                    TT[tmpcc] = PNT_this;
                    //printf("\t%f, %f, %f\n", PNT_this.x, PNT_this.y, PNT_this.z);
                    tmpcc++;
                }
            }
            //printf("------\n\n");
        }
    }

    verts->NumVertsTruncated = tmpcc;

    //printf("[");
    for (int i = 0; i < tmpcc; ++i)
    {
        verts->Verts3DTruncated[i] = TT[i];
        //printf("%f, %f, %f;\n", TT[i].x, TT[i].y, TT[i].z);
    };
    //printf("]");
    return if_touch;
}; // TruncateFracture
template __device__ __host__ bool cuDFNsys::TruncateFracture<double>(cuDFNsys::Fracture<double> *verts,
                                                                     cuDFNsys::Vector1<double> L,
                                                                     int plane,
                                                                     int dir);
template __device__ __host__ bool cuDFNsys::TruncateFracture<float>(cuDFNsys::Fracture<float> *verts,
                                                                    cuDFNsys::Vector1<float> L,
                                                                    int plane,
                                                                    int dir);