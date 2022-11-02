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

#include "Fractures/IdentifyIntersectionKernel.cuh"
// ====================================================
// NAME:        IdentifyIntersectionKernel
// DESCRIPTION: Identify intersection on the GPU side
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================
template <typename T>
__global__ void cuDFNsys::IdentifyIntersectionKernel(cuDFNsys::Fracture<T> *verts,
                                                     int count,
                                                     cuDFNsys::Intersection<T> *Int_sec,
                                                     bool if_trucncated)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;

    if (idx > count - 1)
        return;

    int i = Int_sec[idx].FracIDPair.x;
    int j = Int_sec[idx].FracIDPair.y;
    //Int_sec[idx].IfIntersect = false; // default is false

    T Mat_2_to_3[3][3];
    verts[i].RoationMatrix(Mat_2_to_3, 23); //= verts[i].Roation_Matrix_2Dto3D;
    T Mat_3_to_2[3][3];
    verts[i].RoationMatrix(Mat_3_to_2, 32);

    int NUM_verts = 0;

    if (if_trucncated == false)
        NUM_verts = 4;
    else
        NUM_verts = verts[i].NumVertsTruncated;

    cuDFNsys::Vector3<T> *Frac_verts = (cuDFNsys::Vector3<T> *)malloc(NUM_verts * sizeof(cuDFNsys::Vector3<T>));
    // = verts[i].Verts3D;

    if (if_trucncated == false)
        memcpy(Frac_verts, verts[i].Verts3D, sizeof(cuDFNsys::Vector3<T>) * NUM_verts);
    else
        memcpy(Frac_verts, verts[i].Verts3DTruncated, sizeof(cuDFNsys::Vector3<T>) * NUM_verts);

    cuDFNsys::Vector3<T> Center_ = verts[i].Center;
    cuDFNsys::Vector2<T> *Frac_verts_2D = (cuDFNsys::Vector2<T> *)malloc(NUM_verts * sizeof(cuDFNsys::Vector2<T>));

    for (int k = 0; k < NUM_verts; ++k)
    {
        Frac_verts[k] = cuDFNsys::MakeVector3<T>(Frac_verts[k].x - Center_.x,
                                                 Frac_verts[k].y - Center_.y,
                                                 Frac_verts[k].z - Center_.z);
        cuDFNsys::Vector3<T> HJ = cuDFNsys::ProductSquare3Vector3<T>(Mat_3_to_2, Frac_verts[k]);
        Frac_verts_2D[k] = cuDFNsys::MakeVector2<T>(HJ.x, HJ.y);
    }

    int NUM_verts_j = 0;
    if (if_trucncated == false)
        NUM_verts_j = 4;
    else
        NUM_verts_j = verts[j].NumVertsTruncated;

    cuDFNsys::Vector3<T> *Frac_verts_j = (cuDFNsys::Vector3<T> *)malloc(sizeof(cuDFNsys::Vector3<T>) * NUM_verts_j); // = verts[j].verts_3D;

    if (if_trucncated == false)
        memcpy(Frac_verts_j, verts[j].Verts3D, sizeof(cuDFNsys::Vector3<T>) * NUM_verts_j);
    else
        memcpy(Frac_verts_j, verts[j].Verts3DTruncated, sizeof(cuDFNsys::Vector3<T>) * NUM_verts_j);

    for (int k = 0; k < NUM_verts_j; ++k)
    {
        Frac_verts_j[k] = cuDFNsys::MakeVector3<T>(Frac_verts_j[k].x - Center_.x,
                                                   Frac_verts_j[k].y - Center_.y,
                                                   Frac_verts_j[k].z - Center_.z);
        Frac_verts_j[k] = cuDFNsys::ProductSquare3Vector3<T>(Mat_3_to_2,
                                                             Frac_verts_j[k]);
    }

    //---------------start to identify
    cuDFNsys::Vector3<T> Intersection_xyPlane_Poly2[2];
    bool ik = cuDFNsys::Intersection3DPolyXYPlane<T>(Frac_verts_j,
                                                     NUM_verts_j,
                                                     Intersection_xyPlane_Poly2,
                                                     _TOL_Intersection3DPolyXYPlane);
    if (ik == false)
    {
        free(Frac_verts);
        free(Frac_verts_j);
        free(Frac_verts_2D);
        Int_sec[idx].FracIDPair.x = -1;
        return;
    }
    cuDFNsys::Vector2<T> LineInf[2];
    LineInf[0] = cuDFNsys::MakeVector2<T>(Intersection_xyPlane_Poly2[0].x, Intersection_xyPlane_Poly2[0].y);
    LineInf[1] = cuDFNsys::MakeVector2<T>(Intersection_xyPlane_Poly2[1].x, Intersection_xyPlane_Poly2[1].y);

    cuDFNsys::Vector2<T> Intersection_Poly1_InfLine[2];
    bool jt = cuDFNsys::Intersection2DLine2DPoly<T>(Frac_verts_2D,
                                                    NUM_verts,
                                                    LineInf,
                                                    Intersection_Poly1_InfLine,
                                                    _TOL_Intersection2DLine2DPoly);
    if (jt == false)
    {
        free(Frac_verts);
        free(Frac_verts_j);
        free(Frac_verts_2D);
        Int_sec[idx].FracIDPair.x = -1;
        return;
    }
    int sign_ = -1;
    cuDFNsys::Vector2<T> Intersection_two_segs[2];
    bool kh = cuDFNsys::IntersectionTwoCollinear2DSegs<T>(LineInf,
                                                          Intersection_Poly1_InfLine,
                                                          Intersection_two_segs,
                                                          &sign_,
                                                          _TOL_IntersectionTwoCollinearSegs);

    if (kh == true)
    {
        if (sign_ == 2)
        {
            cuDFNsys::Vector3<T> Intersection_f[2];
            Intersection_f[0] = cuDFNsys::MakeVector3<T>(Intersection_two_segs[0].x,
                                                         Intersection_two_segs[0].y, (T)0);
            Intersection_f[1] = cuDFNsys::MakeVector3<T>(Intersection_two_segs[1].x,
                                                         Intersection_two_segs[1].y, (T)0);

            for (int k = 0; k < 2; ++k)
            {
                Intersection_f[k] = cuDFNsys::ProductSquare3Vector3<T>(Mat_2_to_3,
                                                                       Intersection_f[k]);
                Intersection_f[k] = cuDFNsys::MakeVector3<T>(Intersection_f[k].x + Center_.x, Intersection_f[k].y + Center_.y, Intersection_f[k].z + Center_.z);
            }

            // pair<size_t, size_t> key_ = std::make_pair((size_t)i, (size_t)j);
            // pair<cuDFNsys::Vector3<T>, cuDFNsys::Vector3<T>> p;
            // p.first = Intersection_f[0];
            // p.second = Intersection_f[1];
            // pair<pair<size_t, size_t>, pair<cuDFNsys::Vector3<T>, cuDFNsys::Vector3<T>>> element_ =
            //     std::make_pair(key_, p);
            // Intersection_map.insert(element_);
            Int_sec[idx].Coord[0] = Intersection_f[0];
            Int_sec[idx].Coord[1] = Intersection_f[1];
            free(Frac_verts);
            free(Frac_verts_j);
            free(Frac_verts_2D);
            return;
        }
        else
        {

            free(Frac_verts);
            free(Frac_verts_j);
            free(Frac_verts_2D);
            Int_sec[idx].FracIDPair.x = -1;
            return;
        }
    }
    else
    {

        free(Frac_verts);
        free(Frac_verts_j);
        free(Frac_verts_2D);
        Int_sec[idx].FracIDPair.x = -1;
        return;
    }
}; // IdentifyIntersectionKernel
template __global__ void cuDFNsys::IdentifyIntersectionKernel<double>(cuDFNsys::Fracture<double> *verts,
                                                                      int count,
                                                                      cuDFNsys::Intersection<double> *Int_sec,
                                                                      bool if_trucncated);
template __global__ void cuDFNsys::IdentifyIntersectionKernel<float>(cuDFNsys::Fracture<float> *verts,
                                                                     int count,
                                                                     cuDFNsys::Intersection<float> *Int_sec,
                                                                     bool if_trucncated);