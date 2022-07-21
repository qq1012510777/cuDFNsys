///////////////////////////////////////////////////////////////////
// NAME:              Fracture.cuh
//
// PURPOSE:           A struct of fracture attributes
//
// FUNCTIONS/OBJECTS: Fracture
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../DataTypeSelector/DataTypeSelector.cuh"
#include "../GlobalDef/GlobalDef.cuh"
#include "../Quaternion/Quaternion.cuh"

namespace cuDFNsys
{
template <typename T>
struct Fracture
{
public:
    // conductivity of the fracture
    cuDFNsys::Vector1<T> Conductivity;
    // 3D verts
    cuDFNsys::Vector3<T> Verts3D[4];
    // center
    cuDFNsys::Vector3<T> Center;
    // verts of 3D square/fracture
    cuDFNsys::Vector3<T> Verts3DTruncated[8];
    // number of verts of truncated fractures
    int NumVertsTruncated;
    // radius of circumscribe circle
    cuDFNsys::Vector1<T> Radius;
    // if the fracture intersects the six faces of cubic model?
    bool ConnectModelSurf[6]; // x_min, x_max, y_min, y_max, z_min, z_max,
    // normal vec (normalized)
    cuDFNsys::Vector3<T> NormalVec;

public:
    // get theta value
    __device__ __host__ cuDFNsys::Vector1<T> Theta();
    // get phi value
    __device__ __host__ cuDFNsys::Vector1<T> Phi();
    // get rotation matrix from 2(3) to 3(2)
    __device__ __host__ void RoationMatrix(cuDFNsys::Vector1<T> tmp_R_1[3][3], const int mode);
    // generate 2D verts
    __device__ __host__ void Generate2DVerts(cuDFNsys::Vector2<T> *verts2DDD, uint NUM_verts, bool IfTrimed);
};
}; // namespace cuDFNsys