#include "Geometry/3D/DistancePnt3DPlane.cuh"

// ====================================================
// NAME:        DistancePnt3DPlane
// DESCRIPTION: Distance between a 3D point and a 3D plane
// AUTHOR:      Tingchang YIN
// DATE:        08/04/2022
// ====================================================
__device__ __host__ float cuDFNsys::DistancePnt3DPlane(float3 Plane[3], float3 pnt)
{
    float a, b, c, d;

    float x1 = Plane[0].x;
    float x2 = Plane[1].x;
    float x3 = Plane[2].x;

    float y1 = Plane[0].y;
    float y2 = Plane[1].y;
    float y3 = Plane[2].y;

    float z1 = Plane[0].z;
    float z2 = Plane[1].z;
    float z3 = Plane[2].z;

    float a1 = x2 - x1;
    float b1 = y2 - y1;
    float c1 = z2 - z1;
    float a2 = x3 - x1;
    float b2 = y3 - y1;
    float c2 = z3 - z1;

    a = b1 * c2 - b2 * c1;
    b = a2 * c1 - a1 * c2;
    c = a1 * b2 - b1 * a2;
    d = (-a * x1 - b * y1 - c * z1);

    d = fabs((a * pnt.x + b * pnt.y +
              c * pnt.z + d));

    float e = sqrt(a * a + b * b + c * c);

    return (d / e);
}; // DistancePnt3DPlane