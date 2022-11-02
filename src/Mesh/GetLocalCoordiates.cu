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

#include "Mesh/GetLocalCoordiates.cuh"

// ====================================================
// NAME:        GetLocalCoordiates
// DESCRIPTION: Get 2D local coordiates of elements
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================
template <typename T>
__global__ void cuDFNsys::GetLocalCoordiates(uint3 *element_3D_dev_ptr,
                                             cuDFNsys::Fracture<T> *Frac_verts_device_ptr,
                                             uint *element_Frac_Tag_dev_ptr,
                                             cuDFNsys::EleCoor<T> *coordinate_2D_dev_ptr,
                                             cuDFNsys::Vector3<T> *coordinate_3D_dev_ptr,
                                             int ele_count)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;

    if (i > ele_count - 1)
        return;

    uint FracTag = element_Frac_Tag_dev_ptr[i];

    cuDFNsys::Vector3<T> Center_ = Frac_verts_device_ptr[FracTag].Center;

    T MAT_3to2[3][3];
    Frac_verts_device_ptr[FracTag].RoationMatrix(MAT_3to2, 32);
    //memcpy(MAT_3to2, Frac_verts_device_ptr[FracTag].Roation_Matrix_3Dto2D, sizeof(float) * 9);

    uint node1 = element_3D_dev_ptr[i].x;
    uint node2 = element_3D_dev_ptr[i].y;
    uint node3 = element_3D_dev_ptr[i].z;

    cuDFNsys::Vector3<T> grid_verts[3];
    grid_verts[0] = coordinate_3D_dev_ptr[node1 - 1];
    grid_verts[1] = coordinate_3D_dev_ptr[node2 - 1];
    grid_verts[2] = coordinate_3D_dev_ptr[node3 - 1];

    cuDFNsys::Vector2<T> verts2Dlocal[3];
    for (int j = 0; j < 3; ++j)
    {
        grid_verts[j] = cuDFNsys::MakeVector3(grid_verts[j].x - Center_.x,
                                              grid_verts[j].y - Center_.y,
                                              grid_verts[j].z - Center_.z);

        grid_verts[j] = cuDFNsys::ProductSquare3Vector3<T>(MAT_3to2, grid_verts[j]);

        //coordinate_2D_dev_ptr[i].x[j] = grid_verts[j].x;
        //coordinate_2D_dev_ptr[i].y[j] = grid_verts[j].y;
        verts2Dlocal[j].x = grid_verts[j].x;
        verts2Dlocal[j].y = grid_verts[j].y;
    }

    //coordinate_2D_dev_ptr[i].x[0] = verts2Dlocal[0].x;
    //coordinate_2D_dev_ptr[i].y[0] = verts2Dlocal[0].y;
    //coordinate_2D_dev_ptr[i].x[1] = verts2Dlocal[1].x;
    //coordinate_2D_dev_ptr[i].y[1] = verts2Dlocal[1].y;
    //coordinate_2D_dev_ptr[i].x[2] = verts2Dlocal[2].x;
    //coordinate_2D_dev_ptr[i].y[2] = verts2Dlocal[2].y;

    //-----------check if the triangle orientation with local coordinates is counterclockwise

    bool ori = cuDFNsys::Triangle2DOrientation<T>(verts2Dlocal[0],
                                                  verts2Dlocal[1],
                                                  verts2Dlocal[2]);

    if (ori == true)
    {
        element_3D_dev_ptr[i].x = node1;
        element_3D_dev_ptr[i].y = node3;
        element_3D_dev_ptr[i].z = node2;

        coordinate_2D_dev_ptr[i].x[0] = verts2Dlocal[0].x;
        coordinate_2D_dev_ptr[i].y[0] = verts2Dlocal[0].y;

        coordinate_2D_dev_ptr[i].x[1] = verts2Dlocal[2].x;
        coordinate_2D_dev_ptr[i].y[1] = verts2Dlocal[2].y;

        coordinate_2D_dev_ptr[i].x[2] = verts2Dlocal[1].x;
        coordinate_2D_dev_ptr[i].y[2] = verts2Dlocal[1].y;

        //printf("clockwise: %d %d %d, change to: %d %d %d\ncoord: (%lf %lf), (%lf %lf), (%lf %lf) changed to (%lf %lf), (%lf %lf), (%lf %lf)\n",
        //       node1, node2, node3, element_3D_dev_ptr[i].x, element_3D_dev_ptr[i].y, element_3D_dev_ptr[i].z,
        //       verts2Dlocal[0].x, verts2Dlocal[0].y, verts2Dlocal[1].x, verts2Dlocal[1].y, verts2Dlocal[2].x, verts2Dlocal[2].y,
        //       coordinate_2D_dev_ptr[i].x[0],
        //       coordinate_2D_dev_ptr[i].y[0],
        //       coordinate_2D_dev_ptr[i].x[1],
        //       coordinate_2D_dev_ptr[i].y[1],
        //       coordinate_2D_dev_ptr[i].x[2],
        //       coordinate_2D_dev_ptr[i].y[2]);
    }
    else
    {
        element_3D_dev_ptr[i].x = node1;
        element_3D_dev_ptr[i].y = node2;
        element_3D_dev_ptr[i].z = node3;

        coordinate_2D_dev_ptr[i].x[0] = verts2Dlocal[0].x;
        coordinate_2D_dev_ptr[i].y[0] = verts2Dlocal[0].y;

        coordinate_2D_dev_ptr[i].x[1] = verts2Dlocal[1].x;
        coordinate_2D_dev_ptr[i].y[1] = verts2Dlocal[1].y;

        coordinate_2D_dev_ptr[i].x[2] = verts2Dlocal[2].x;
        coordinate_2D_dev_ptr[i].y[2] = verts2Dlocal[2].y;
    }
}; //GetLocalCoordiates
template __global__ void cuDFNsys::GetLocalCoordiates<double>(uint3 *element_3D_dev_ptr,
                                                              cuDFNsys::Fracture<double> *Frac_verts_device_ptr,
                                                              uint *element_Frac_Tag_dev_ptr,
                                                              cuDFNsys::EleCoor<double> *coordinate_2D_dev_ptr,
                                                              cuDFNsys::Vector3<double> *coordinate_3D_dev_ptr,
                                                              int ele_count);
template __global__ void cuDFNsys::GetLocalCoordiates<float>(uint3 *element_3D_dev_ptr,
                                                             cuDFNsys::Fracture<float> *Frac_verts_device_ptr,
                                                             uint *element_Frac_Tag_dev_ptr,
                                                             cuDFNsys::EleCoor<float> *coordinate_2D_dev_ptr,
                                                             cuDFNsys::Vector3<float> *coordinate_3D_dev_ptr,
                                                             int ele_count);