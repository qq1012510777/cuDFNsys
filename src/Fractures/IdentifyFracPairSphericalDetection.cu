#include "Fractures/IdentifyFracPairSphericalDetection.cuh"

// ====================================================
// NAME:        IdentifyFracPairSphericalDetection
// DESCRIPTION: Identify fracture pair where the two circumscribed
//              spheres intersect
// AUTHOR:      Tingchang YIN
// DATE:        09/04/2022
// ====================================================
__global__ void cuDFNsys::IdentifyFracPairSphericalDetection(cuDFNsys::Fracture *verts,
                                                             int3 *Frac_pairs,
                                                             int InitialPairNO,
                                                             int count)
{
    int idx_TT = threadIdx.x + blockIdx.x * blockDim.x;

    if (idx_TT > count - 1)
        return;

    int idx = idx_TT + InitialPairNO;

    int x_ = floor((pow(2 * (idx + 1), 0.5) + 1 / 2.0));
    int y_ = idx - 0.5 * x_ * (x_ - 1);
    //printf("%d: x_ = %d, y_ =  %d\n",idx, x_, y_);
    Frac_pairs[idx_TT].x = x_;
    Frac_pairs[idx_TT].y = y_;
    Frac_pairs[idx_TT].z = 1;

    float3 dist_two_frac = make_float3(verts[x_].Center.x - verts[y_].Center.x,
                                       verts[x_].Center.y - verts[y_].Center.y,
                                       verts[x_].Center.z - verts[y_].Center.z);

    float ddis = pow(dist_two_frac.x * dist_two_frac.x +
                         dist_two_frac.y * dist_two_frac.y +
                         dist_two_frac.z * dist_two_frac.z,
                     0.5);

    if (ddis > (verts[x_].Radius + verts[y_].Radius))
        Frac_pairs[idx_TT].z = 0;
};