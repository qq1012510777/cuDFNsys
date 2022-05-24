#include "Mesh/IdentifyEleFrac.cuh"

// ====================================================
// NAME:        IdentifyEleFrac
// DESCRIPTION: Identify fracture ID of elements 
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================
__global__ void cuDFNsys::IdentifyEleFrac(uint3 *One_entity_one_ele_dev_ptr,
                                          float3 *coordinate_3D_dev_ptr,
                                          cuDFNsys::Fracture *Frac_verts_device_ptr,
                                          int *Elements_Frac_dev_ptr,
                                          int entity_count,
                                          int frac_count,
                                          float _tol_)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;

    if (i > entity_count - 1)
        return;

    //printf("entering Identify %d\n", i);

    Elements_Frac_dev_ptr[i] = -1;
    //printf("before identify %d\n", Elements_Frac_dev_ptr[i]);

    uint node1 = One_entity_one_ele_dev_ptr[i].x;
    uint node2 = One_entity_one_ele_dev_ptr[i].y;
    uint node3 = One_entity_one_ele_dev_ptr[i].z;
    //printf("node %d %d %d\n", node1, node2, node3);

    float3 vert1 = coordinate_3D_dev_ptr[node1 - 1];
    float3 vert2 = coordinate_3D_dev_ptr[node2 - 1];
    float3 vert3 = coordinate_3D_dev_ptr[node3 - 1];

    //printf("kernel 0, sizeof fracs: %d\n", frac_count);
    for (int j = 0; j < frac_count; ++j)
    {
        //printf("kernel 1\n");

        float3 Plane[3] = {Frac_verts_device_ptr[j].Verts3D[0], Frac_verts_device_ptr[j].Verts3D[1], Frac_verts_device_ptr[j].Verts3D[2]};

        float3 Center_l = (Frac_verts_device_ptr[j].Center);

        float3 *verts_ele[3] = {&(vert1), &(vert2), &(vert3)};

        float MAT_3to2[3][3];
        Frac_verts_device_ptr[j].RoationMatrix(MAT_3to2, 32);
        //memcpy(MAT_3to2, Frac_verts_device_ptr[j].Roation_Matrix_3Dto2D, sizeof(float) * 9);

        float2 verts_2D__s[4];
        Frac_verts_device_ptr[j].Generate2DVerts(verts_2D__s, 4, false);
        //memcpy(verts_2D__s, Frac_verts_device_ptr[j].verts_2D, sizeof(float2) * 4);

        bool belong_to_this_frac = true;
        //printf("kernel 2\n");
        for (int k = 0; k < 3; ++k)
        {
            float dis = cuDFNsys::DistancePnt3DPlane(Plane, (*verts_ele[k]));

            if (abs(dis) > _tol_)
            {
                //printf("%f\n", abs(dis));
                belong_to_this_frac = false;
                break;
            }

            float3 Pnt_l = make_float3((*verts_ele[k]).x - (Center_l).x, (*verts_ele[k]).y - (Center_l).y, (*verts_ele[k]).z - (Center_l).z);
            Pnt_l = ProductSquare3Float3(MAT_3to2, Pnt_l);

            float2 Pnt_ll = make_float2(Pnt_l.x, Pnt_l.y);

            bool inside_ = cuDFNsys::IfPntInside2DConvexPoly(Pnt_ll, verts_2D__s, 4);
            bool on_bound = cuDFNsys::IfPntLiesOnBound2DConvexPoly(Pnt_ll, verts_2D__s, 4, _tol_);
            //printf("inside %d on_bound %d \n", inside_, on_bound);
            if (inside_ == false && on_bound == false)
            {
                belong_to_this_frac = false;

                //printf("entity %d, coplane but not inside! the %d point\n", i, k + 1);
                //printf("\tL_%dpnt(1, :)= [%f %f %f];\n", i, (*verts_ele[0]).x, (*verts_ele[0]).y, (*verts_ele[0]).z);
                //printf("\tL_%dpnt(2, :)= [%f %f %f];\n", i, (*verts_ele[1]).x, (*verts_ele[1]).y, (*verts_ele[1]).z);
                //printf("\tL_%dpnt(3, :)= [%f %f %f];\n", i, (*verts_ele[2]).x, (*verts_ele[2]).y, (*verts_ele[2]).z);
                //
                //printf("\tL_%dverts(1, :)=[%f %f %f];\n", i, Frac_verts_device_ptr[j].Verts3D[0].x, Frac_verts_device_ptr[j].Verts3D[0].y, Frac_verts_device_ptr[j].Verts3D[0].z);
                //printf("\tL_%dverts(2, :)=[%f %f %f];\n", i, Frac_verts_device_ptr[j].Verts3D[1].x, Frac_verts_device_ptr[j].Verts3D[1].y, Frac_verts_device_ptr[j].Verts3D[1].z);
                //printf("\tL_%dverts(3, :)=[%f %f %f];\n", i, Frac_verts_device_ptr[j].Verts3D[2].x, Frac_verts_device_ptr[j].Verts3D[2].y, Frac_verts_device_ptr[j].Verts3D[2].z);
                //printf("\tL_%dverts(4, :)=[%f %f %f];\n", i, Frac_verts_device_ptr[j].Verts3D[3].x, Frac_verts_device_ptr[j].Verts3D[3].y, Frac_verts_device_ptr[j].Verts3D[3].z);
                break;
            }
        }

        if (belong_to_this_frac == true)
        {
            Elements_Frac_dev_ptr[i] = j;
            //printf("entity-%d\n", Elements_Frac_dev_ptr[i]);
            break;
        }
    }

    //printf("after identify %d\n", Elements_Frac_dev_ptr[i]);
}; // IdentifyEleFrac