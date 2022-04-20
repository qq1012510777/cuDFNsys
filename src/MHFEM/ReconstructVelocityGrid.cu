#include "MHFEM/ReconstructVelocityGrid.cuh"

// ====================================================
// NAME:        ReconstructVelocityGrid
// DESCRIPTION: Reconstruct the velocity field in a grid
// AUTHOR:      Tingchang YIN
// DATE:        20/04/2022
// ====================================================
__host__ __device__ float2 cuDFNsys::ReconstructVelocityGrid(float2 Point_,
                                                             float2 Vertex[3],
                                                             float3 VelocityEdgeNormal)
{
    float2 velocity_;

    float2 A = make_float2(Vertex[2].x - Vertex[1].x, Vertex[2].y - Vertex[1].y);
    float2 B = make_float2(Vertex[0].x - Vertex[2].x, Vertex[0].y - Vertex[2].y);
    float2 C = make_float2(Vertex[1].x - Vertex[0].x, Vertex[1].y - Vertex[0].y);

    float3 edge_length = make_float3(sqrt(A.x * A.x + A.y * A.y),
                                     sqrt(B.x * B.x + B.y * B.y),
                                     sqrt(C.x * C.x + C.y * C.y));
    float Area = cuDFNsys::Triangle2DArea(Vertex[0], Vertex[1], Vertex[2]);

    float2 Phi_1 = make_float2(edge_length.x / (2 * Area) * (Point_.x - Vertex[0].x) * VelocityEdgeNormal.y,
                               edge_length.x / (2 * Area) * (Point_.y - Vertex[0].y) * VelocityEdgeNormal.y);

    float2 Phi_2 = make_float2(edge_length.y / (2 * Area) * (Point_.x - Vertex[1].x) * VelocityEdgeNormal.z,
                               edge_length.y / (2 * Area) * (Point_.y - Vertex[1].y) * VelocityEdgeNormal.z);

    float2 Phi_3 = make_float2(edge_length.z / (2 * Area) * (Point_.x - Vertex[2].x) * VelocityEdgeNormal.x,
                               edge_length.z / (2 * Area) * (Point_.y - Vertex[2].y) * VelocityEdgeNormal.x);

    velocity_.x = Phi_1.x + Phi_2.x + Phi_3.x;
    velocity_.y = Phi_1.y + Phi_2.y + Phi_3.y;
    return velocity_;
}; // ReconstructVelocityGrid