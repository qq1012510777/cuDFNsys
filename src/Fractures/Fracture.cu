#include "Fractures/Fracture.cuh"

// ====================================================
// NAME:        Theta
// DESCRIPTION: get theta value (radian).
// AUTHOR:      Tingchang YIN
// DATE:        02/04/2022
// ====================================================
template <typename T>
__device__ __host__ cuDFNsys::Vector1<T> cuDFNsys::Fracture<T>::Theta()
{
    return acos(NormalVec.z);
}; // Theta
template __device__ __host__ cuDFNsys::Vector1<double> cuDFNsys::Fracture<double>::Theta();
template __device__ __host__ cuDFNsys::Vector1<float> cuDFNsys::Fracture<float>::Theta();

// ====================================================
// NAME:        Phi
// DESCRIPTION: get Phi Phi (radian).
// AUTHOR:      Tingchang YIN
// DATE:        02/04/2022
// ====================================================
template <typename T>
__device__ __host__ cuDFNsys::Vector1<T> cuDFNsys::Fracture<T>::Phi()
{
    cuDFNsys::Vector1<T> phi = atan2(NormalVec.y, NormalVec.x);
    return (phi > 0 ? phi : phi + 2.0 * M_PI);
}; // Phi
template __device__ __host__ cuDFNsys::Vector1<double> cuDFNsys::Fracture<double>::Phi();
template __device__ __host__ cuDFNsys::Vector1<float> cuDFNsys::Fracture<float>::Phi();

// ====================================================
// NAME:        RoationMatrix
// DESCRIPTION: get RoationMatrix from 3(2) to 2(3).
//              mode is 32 or 23
// AUTHOR:      Tingchang YIN
// DATE:        02/04/2022
// ====================================================
template <typename T>
__device__ __host__ void cuDFNsys::Fracture<T>::RoationMatrix(cuDFNsys::Vector1<T> tmp_R_1[3][3], const int mode)
{
    cuDFNsys::Vector3<T> rotate_axis;
    rotate_axis.x = -NormalVec.y, rotate_axis.y = NormalVec.x, rotate_axis.z = 0.00;

    cuDFNsys::Vector1<T> norm_axis = sqrt(rotate_axis.x * rotate_axis.x + rotate_axis.y * rotate_axis.y);
    rotate_axis.x /= norm_axis;
    rotate_axis.y /= norm_axis;

    int sign_ = 0;
    if (mode == 32)
        sign_ = -1;
    else if (mode == 23)
        sign_ = 1;

    cuDFNsys::Quaternion<T> Qua;
    Qua = Qua.DescribeRotation(rotate_axis, sign_ * acos(NormalVec.z));

    cuDFNsys::Vector4<T> quater_ = Qua.GetQuaternionNum();
    cuDFNsys::Vector1<T> w = quater_.x,
                         x = quater_.y,
                         y = quater_.z,
                         z = quater_.w;
    cuDFNsys::Vector1<T> tmp_R[3][3] = {1 - 2 * y * y - 2 * z * z, 2 * x * y - 2 * w * z, 2 * x * z + 2 * w * y,
                                        2 * x * y + 2 * w * z, 1 - 2 * x * x - 2 * z * z, 2 * y * z - 2 * w * x,
                                        2 * x * z - 2 * w * y, 2 * y * z + 2 * w * x, 1 - 2 * x * x - 2 * y * y};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            tmp_R_1[i][j] = tmp_R[i][j];
}; // RoationMatrix
template __device__ __host__ void cuDFNsys::Fracture<double>::RoationMatrix(cuDFNsys::Vector1<double> tmp_R_1[3][3], const int mode);
template __device__ __host__ void cuDFNsys::Fracture<float>::RoationMatrix(cuDFNsys::Vector1<float> tmp_R_1[3][3], const int mode);

// ====================================================
// NAME:        Generate2DVerts
// DESCRIPTION: generate 2D verts.
// AUTHOR:      Tingchang YIN
// DATE:        02/04/2022
// ====================================================
template <typename T>
__device__ __host__ void cuDFNsys::Fracture<T>::Generate2DVerts(cuDFNsys::Vector2<T> *verts2DDD, uint NUM_verts, bool IfTrimed)
{
    cuDFNsys::Vector3<T> rotate_axis;
    rotate_axis.x = -NormalVec.y, rotate_axis.y = NormalVec.x, rotate_axis.z = 0.00;

    cuDFNsys::Vector1<T> norm_axis = sqrt(rotate_axis.x * rotate_axis.x + rotate_axis.y * rotate_axis.y);
    rotate_axis.x /= norm_axis;
    rotate_axis.y /= norm_axis;
    cuDFNsys::Quaternion<T> Qua;
    Qua = Qua.DescribeRotation(rotate_axis, -1.0 * acos(this->NormalVec.z));

    if (IfTrimed == false)
        NUM_verts = 4;
    for (uint i = 0; i < NUM_verts; ++i)
    {
        cuDFNsys::Vector3<T> Vertex__;

        if (IfTrimed == false)
        {
            Vertex__.x = this->Verts3D[i].x - this->Center.x,
            Vertex__.y = this->Verts3D[i].y - this->Center.y,
            Vertex__.z = this->Verts3D[i].z - this->Center.z;
        }
        else
        {
            Vertex__.x = this->Verts3DTruncated[i].x - this->Center.x,
            Vertex__.y = this->Verts3DTruncated[i].y - this->Center.y,
            Vertex__.z = this->Verts3DTruncated[i].z - this->Center.z;
        }

        Vertex__ = Qua.Rotate(Vertex__);
        verts2DDD[i].x = Vertex__.x;
        verts2DDD[i].y = Vertex__.y;
    };
}; // Generate2DVerts
template __device__ __host__ void cuDFNsys::Fracture<double>::Generate2DVerts(cuDFNsys::Vector2<double> *verts2DDD, uint NUM_verts, bool IfTrimed);
template __device__ __host__ void cuDFNsys::Fracture<float>::Generate2DVerts(cuDFNsys::Vector2<float> *verts2DDD, uint NUM_verts, bool IfTrimed);