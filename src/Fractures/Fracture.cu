#include "Fractures/Fracture.cuh"

// ====================================================
// NAME:        Theta
// DESCRIPTION: get theta value (radian).
// AUTHOR:      Tingchang YIN
// DATE:        02/04/2022
// ====================================================
__device__ __host__ float cuDFNsys::Fracture::Theta()
{
    return acos(NormalVec.z);
}; // Theta

// ====================================================
// NAME:        Phi
// DESCRIPTION: get Phi Phi (radian).
// AUTHOR:      Tingchang YIN
// DATE:        02/04/2022
// ====================================================
__device__ __host__ float cuDFNsys::Fracture::Phi()
{
    float phi = atan2(NormalVec.y, NormalVec.x);
    return (phi > 0 ? phi : phi + 2.0 * M_PI);
}; // Phi

// ====================================================
// NAME:        RoationMatrix
// DESCRIPTION: get RoationMatrix from 3(2) to 2(3).
//              mode is 32 or 23
// AUTHOR:      Tingchang YIN
// DATE:        02/04/2022
// ====================================================
__device__ __host__ void cuDFNsys::Fracture::RoationMatrix(float tmp_R_1[3][3], const int mode)
{
    float3 rotate_axis = make_float3(-NormalVec.y, NormalVec.x, 0.00);
    float norm_axis = sqrt(rotate_axis.x * rotate_axis.x + rotate_axis.y * rotate_axis.y);
    rotate_axis.x /= norm_axis;
    rotate_axis.y /= norm_axis;

    int sign_ = 0;
    if (mode == 32)
        sign_ = -1;
    else if (mode == 23)
        sign_ = 1;

    cuDFNsys::Quaternion Qua;
    Qua = Qua.DescribeRotation(rotate_axis, sign_ * acos(NormalVec.z));

    float4 quater_ = Qua.GetQuaternionNum();
    float w = quater_.x;
    float x = quater_.y;
    float y = quater_.z;
    float z = quater_.w;
    float tmp_R[3][3] = {1 - 2 * y * y - 2 * z * z, 2 * x * y - 2 * w * z, 2 * x * z + 2 * w * y,
                         2 * x * y + 2 * w * z, 1 - 2 * x * x - 2 * z * z, 2 * y * z - 2 * w * x,
                         2 * x * z - 2 * w * y, 2 * y * z + 2 * w * x, 1 - 2 * x * x - 2 * y * y};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            tmp_R_1[i][j] = tmp_R[i][j];
}; // RoationMatrix

// ====================================================
// NAME:        Generate2DVerts
// DESCRIPTION: generate 2D verts.
// AUTHOR:      Tingchang YIN
// DATE:        02/04/2022
// ====================================================
__device__ __host__ void cuDFNsys::Fracture::Generate2DVerts(float2 verts2DDD[4])
{
    float3 rotate_axis = make_float3(-NormalVec.y, NormalVec.x, 0.00);
    float norm_axis = sqrt(rotate_axis.x * rotate_axis.x + rotate_axis.y * rotate_axis.y);
    rotate_axis.x /= norm_axis;
    rotate_axis.y /= norm_axis;
    cuDFNsys::Quaternion Qua;
    Qua = Qua.DescribeRotation(rotate_axis, -1.0 * acos(this->NormalVec.z));

    float2 verts2D[4];

    float3 verts3D__[4];
    for (int i = 0; i < 4; ++i)
    {
        verts3D__[i] = make_float3(this->Verts3D[i].x - this->Center.x,
                                   this->Verts3D[i].y - this->Center.y,
                                   this->Verts3D[i].z - this->Center.z);
        verts3D__[i] = Qua.Rotate(verts3D__[i]);
        verts2D[i].x = verts3D__[i].x;
        verts2D[i].y = verts3D__[i].y;
    }
    for (int i = 0; i < 4; ++i)
        verts2DDD[i] = verts2D[i];
}; // Generate2DVerts