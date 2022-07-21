#include "MHFEM/ReconstructVelocityGrid.cuh"

// ====================================================
// NAME:        ReconstructVelocityGrid
// DESCRIPTION: Reconstruct the velocity field in a grid
// AUTHOR:      Tingchang YIN
// DATE:        20/04/2022
// ====================================================
template <typename T>
__host__ __device__ cuDFNsys::Vector2<T> cuDFNsys::ReconstructVelocityGrid(cuDFNsys::Vector2<T> Point_,
                                                                           cuDFNsys::Vector2<T> Vertex[3],
                                                                           cuDFNsys::Vector3<T> VelocityEdgeNormal)
{
    cuDFNsys::Vector2<T> velocity_;

    cuDFNsys::Vector2<T> A = cuDFNsys::MakeVector2(Vertex[2].x - Vertex[1].x, Vertex[2].y - Vertex[1].y);
    cuDFNsys::Vector2<T> B = cuDFNsys::MakeVector2(Vertex[0].x - Vertex[2].x, Vertex[0].y - Vertex[2].y);
    cuDFNsys::Vector2<T> C = cuDFNsys::MakeVector2(Vertex[1].x - Vertex[0].x, Vertex[1].y - Vertex[0].y);

    cuDFNsys::Vector3<T> edge_length = cuDFNsys::MakeVector3(sqrt(A.x * A.x + A.y * A.y),
                                                             sqrt(B.x * B.x + B.y * B.y),
                                                             sqrt(C.x * C.x + C.y * C.y));
    T Area = cuDFNsys::Triangle2DArea<T>(Vertex[0], Vertex[1], Vertex[2]);

    cuDFNsys::Vector2<T> Phi_1 = cuDFNsys::MakeVector2(edge_length.x / (2 * Area) * (Point_.x - Vertex[0].x) * VelocityEdgeNormal.y,
                                                       edge_length.x / (2 * Area) * (Point_.y - Vertex[0].y) * VelocityEdgeNormal.y);

    cuDFNsys::Vector2<T> Phi_2 = cuDFNsys::MakeVector2(edge_length.y / (2 * Area) * (Point_.x - Vertex[1].x) * VelocityEdgeNormal.z,
                                                       edge_length.y / (2 * Area) * (Point_.y - Vertex[1].y) * VelocityEdgeNormal.z);

    cuDFNsys::Vector2<T> Phi_3 = cuDFNsys::MakeVector2(edge_length.z / (2 * Area) * (Point_.x - Vertex[2].x) * VelocityEdgeNormal.x,
                                                       edge_length.z / (2 * Area) * (Point_.y - Vertex[2].y) * VelocityEdgeNormal.x);

    velocity_.x = Phi_1.x + Phi_2.x + Phi_3.x;
    velocity_.y = Phi_1.y + Phi_2.y + Phi_3.y;
    return velocity_;
}; // ReconstructVelocityGrid
template __host__ __device__ cuDFNsys::Vector2<double> cuDFNsys::ReconstructVelocityGrid<double>(cuDFNsys::Vector2<double> Point_,
                                                                                                 cuDFNsys::Vector2<double> Vertex[3],
                                                                                                 cuDFNsys::Vector3<double> VelocityEdgeNormal);
template __host__ __device__ cuDFNsys::Vector2<float> cuDFNsys::ReconstructVelocityGrid<float>(cuDFNsys::Vector2<float> Point_,
                                                                                               cuDFNsys::Vector2<float> Vertex[3],
                                                                                               cuDFNsys::Vector3<float> VelocityEdgeNormal);