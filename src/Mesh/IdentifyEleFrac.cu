#include "Mesh/IdentifyEleFrac.cuh"

// ====================================================
// NAME:        IdentifyEleFrac
// DESCRIPTION: Identify fracture ID of elements
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================
template <typename T>
__global__ void cuDFNsys::IdentifyEleFrac(uint3 *One_entity_one_ele_dev_ptr,
                                          cuDFNsys::Vector3<T> *coordinate_3D_dev_ptr,
                                          cuDFNsys::Fracture<T> *Frac_verts_device_ptr,
                                          int *Elements_Frac_dev_ptr,
                                          int entity_count,
                                          int frac_count,
                                          T _tol_)
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

    cuDFNsys::Vector3<T> vert1 = coordinate_3D_dev_ptr[node1 - 1];
    cuDFNsys::Vector3<T> vert2 = coordinate_3D_dev_ptr[node2 - 1];
    cuDFNsys::Vector3<T> vert3 = coordinate_3D_dev_ptr[node3 - 1];

    //printf("kernel 0, sizeof fracs: %d\n", frac_count);
    for (int j = 0; j < frac_count; ++j)
    {
        //printf("kernel 1\n");

        cuDFNsys::Vector3<T> Plane[3] = {Frac_verts_device_ptr[j].Verts3D[0],
                                         Frac_verts_device_ptr[j].Verts3D[1],
                                         Frac_verts_device_ptr[j].Verts3D[2]};

        cuDFNsys::Vector3<T> Center_l = (Frac_verts_device_ptr[j].Center);

        cuDFNsys::Vector3<T> *verts_ele[3] = {&(vert1), &(vert2), &(vert3)};

        T MAT_3to2[3][3];
        Frac_verts_device_ptr[j].RoationMatrix(MAT_3to2, 32);
        //memcpy(MAT_3to2, Frac_verts_device_ptr[j].Roation_Matrix_3Dto2D, sizeof(float) * 9);

        cuDFNsys::Vector2<T> verts_2D__s[4];
        Frac_verts_device_ptr[j].Generate2DVerts(verts_2D__s, 4, false);
        //memcpy(verts_2D__s, Frac_verts_device_ptr[j].verts_2D, sizeof(float2) * 4);

        bool belong_to_this_frac = true;
        //printf("kernel 2\n");
        for (int k = 0; k < 3; ++k)
        {
            T dis = cuDFNsys::DistancePnt3DPlane<T>(Plane, (*verts_ele[k]));

            if (abs(dis) > _tol_)
            {
                //printf("%f\n", abs(dis));
                belong_to_this_frac = false;
                break;
            }

            cuDFNsys::Vector3<T> Pnt_l = cuDFNsys::MakeVector3((*verts_ele[k]).x - (Center_l).x,
                                                               (*verts_ele[k]).y - (Center_l).y,
                                                               (*verts_ele[k]).z - (Center_l).z);
            Pnt_l = cuDFNsys::ProductSquare3Vector3<T>(MAT_3to2, Pnt_l);

            cuDFNsys::Vector2<T> Pnt_ll = cuDFNsys::MakeVector2(Pnt_l.x, Pnt_l.y);

            bool inside_ = cuDFNsys::IfPntInside2DConvexPoly<T>(Pnt_ll, verts_2D__s, 4);
            bool on_bound = cuDFNsys::IfPntLiesOnBound2DConvexPoly<T>(Pnt_ll, verts_2D__s, 4, _tol_);
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
template __global__ void cuDFNsys::IdentifyEleFrac<double>(uint3 *One_entity_one_ele_dev_ptr,
                                                           cuDFNsys::Vector3<double> *coordinate_3D_dev_ptr,
                                                           cuDFNsys::Fracture<double> *Frac_verts_device_ptr,
                                                           int *Elements_Frac_dev_ptr,
                                                           int entity_count,
                                                           int frac_count,
                                                           double _tol_);
template __global__ void cuDFNsys::IdentifyEleFrac<float>(uint3 *One_entity_one_ele_dev_ptr,
                                                          cuDFNsys::Vector3<float> *coordinate_3D_dev_ptr,
                                                          cuDFNsys::Fracture<float> *Frac_verts_device_ptr,
                                                          int *Elements_Frac_dev_ptr,
                                                          int entity_count,
                                                          int frac_count,
                                                          float _tol_);