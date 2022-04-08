#include "Mesh/GetLocalCoordiates.cuh"

// ====================================================
// NAME:        GetLocalCoordiates
// DESCRIPTION: Get 2D local coordiates of elements 
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================
__global__ void cuDFNsys::GetLocalCoordiates(uint3 *element_3D_dev_ptr,
                                             cuDFNsys::Fracture *Frac_verts_device_ptr,
                                             uint *element_Frac_Tag_dev_ptr,
                                             cuDFNsys::EleCoor *coordinate_2D_dev_ptr,
                                             float3 *coordinate_3D_dev_ptr,
                                             int ele_count)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;

    if (i > ele_count - 1)
        return;

    uint FracTag = element_Frac_Tag_dev_ptr[i];

    float3 Center_ = Frac_verts_device_ptr[FracTag].Center;

    float MAT_3to2[3][3];
    Frac_verts_device_ptr[FracTag].RoationMatrix(MAT_3to2, 32);
    //memcpy(MAT_3to2, Frac_verts_device_ptr[FracTag].Roation_Matrix_3Dto2D, sizeof(float) * 9);

    uint node1 = element_3D_dev_ptr[i].x;
    uint node2 = element_3D_dev_ptr[i].y;
    uint node3 = element_3D_dev_ptr[i].z;

    float3 grid_verts[3];
    grid_verts[0] = coordinate_3D_dev_ptr[node1 - 1];
    grid_verts[1] = coordinate_3D_dev_ptr[node2 - 1];
    grid_verts[2] = coordinate_3D_dev_ptr[node3 - 1];

    float2 verts2Dlocal[3];
    for (int j = 0; j < 3; ++j)
    {
        grid_verts[j] = make_float3(grid_verts[j].x - Center_.x,
                                    grid_verts[j].y - Center_.y,
                                    grid_verts[j].z - Center_.z);

        grid_verts[j] = cuDFNsys::ProductSquare3Float3(MAT_3to2, grid_verts[j]);

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

    bool ori = cuDFNsys::Triangle2DOrientation(verts2Dlocal[0],
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