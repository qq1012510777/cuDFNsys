#include "Geometry/2D/IfPntInside2DConvexPoly.cuh"
// ====================================================
// NAME:        IfPntInside2DConvexPoly
// DESCRIPTION: if a point is inside a 2D polygon
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================
__device__ __host__ bool cuDFNsys::IfPntInside2DConvexPoly(float2 pnt, float2 *verts, int N)
{
    bool inside = true;

    float R[10];

    for (int i = 0; i < N; ++i)
    {
        float2 A = make_float2(verts[i].x - pnt.x, verts[i].y - pnt.y);
        float norm_A = pow(A.x * A.x + A.y * A.y, 0.5);
        A = make_float2(A.x / norm_A, A.y / norm_A);

        float2 B = make_float2(verts[(i + 1) % N].x - pnt.x, verts[(i + 1) % N].y - pnt.y);
        float norm_B = pow(B.x * B.x + B.y * B.y, 0.5);
        B = make_float2(B.x / norm_B, B.y / norm_B);

        R[i] = cuDFNsys::CrossProductFloat2(A, B);

        R[i] = R[i] / abs(R[i]);
    };

    for (int i = 0; i < N; ++i)
    {
        if (R[i] != R[(i + 1) % N])
        {
            inside = false;
            break;
        }
    }

    return inside;
}; // IfPntInside2DConvexPoly