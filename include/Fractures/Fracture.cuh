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
#include "../GlobalDef/GlobalDef.cuh"
#include "../Quaternion/Quaternion.cuh"

namespace cuDFNsys
{
struct Fracture
{
public:
    // conductivity of the fracture
    float Conductivity;
    // 3D verts
    float3 Verts3D[4];
    // center
    float3 Center;
    // verts of 3D square/fracture
    float3 Verts3DTruncated[8];
    // number of verts of truncated fractures
    int NumVertsTruncated;
    // radius of circumscribe circle
    float Radius;
    // if the fracture intersects the six faces of cubic model?
    bool ConnectModelSurf[6]; // x_min, x_max, y_min, y_max, z_min, z_max,
    // normal vec (normalized)
    float3 NormalVec;

public:
    // get theta value
    __device__ __host__ float Theta();
    // get phi value
    __device__ __host__ float Phi();
    // get rotation matrix from 2(3) to 3(2)
    __device__ __host__ void RoationMatrix(float tmp_R_1[3][3], const int mode);
    // generate 2D verts
    __device__ __host__ void Generate2DVerts(float2 verts2DDD[4]);
};
}; // namespace cuDFNsys