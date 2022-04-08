#include "Fractures/IdentifyIntersection.cuh"

// ====================================================
// NAME:        IdentifyIntersection
// DESCRIPTION: Identify Intersections of fractures
//              in a  DFN with CPU
// AUTHOR:      Tingchang YIN
// DATE:        07/04/2022
// ====================================================
cuDFNsys::IdentifyIntersection::IdentifyIntersection(thrust::host_vector<cuDFNsys::Fracture> verts,
                                                     const bool &if_trucncated,
                                                     MapIntersection &Intersection_map)
{
    Intersection_map.clear();

    for (int i = 1; i < verts.size(); ++i)
    {
        for (int j = 0; j < i; ++j)
        {
            float3 dist_two_frac = make_float3(verts[i].Center.x - verts[j].Center.x,
                                               verts[i].Center.y - verts[j].Center.y,
                                               verts[i].Center.z - verts[j].Center.z);

            float ddis = pow(dist_two_frac.x * dist_two_frac.x +
                                 dist_two_frac.y * dist_two_frac.y +
                                 dist_two_frac.z * dist_two_frac.z,
                             0.5);

            if (ddis > (verts[i].Radius + verts[j].Radius))
                continue;
            //---------------

            float Mat_2_to_3[3][3];
            verts[i].RoationMatrix(Mat_2_to_3, 23); //= verts[i].Roation_Matrix_2Dto3D;
            float Mat_3_to_2[3][3];
            verts[i].RoationMatrix(Mat_3_to_2, 32);

            int NUM_verts = 0;

            if (if_trucncated == false)
                NUM_verts = 4;
            else
                NUM_verts = verts[i].NumVertsTruncated;

            float3 *Frac_verts = (float3 *)malloc(NUM_verts * sizeof(float3));
            // = verts[i].Verts3D;

            if (if_trucncated == false)
                memcpy(Frac_verts, verts[i].Verts3D, sizeof(float3) * NUM_verts);
            else
                memcpy(Frac_verts, verts[i].Verts3DTruncated, sizeof(float3) * NUM_verts);

            float3 Center_ = verts[i].Center;
            float2 *Frac_verts_2D = (float2 *)malloc(NUM_verts * sizeof(float2));

            for (int k = 0; k < NUM_verts; ++k)
            {
                Frac_verts[k] = make_float3(Frac_verts[k].x - Center_.x,
                                            Frac_verts[k].y - Center_.y,
                                            Frac_verts[k].z - Center_.z);
                float3 HJ = cuDFNsys::ProductSquare3Float3(Mat_3_to_2, Frac_verts[k]);
                Frac_verts_2D[k] = make_float2(HJ.x, HJ.y);
            }

            int NUM_verts_j = 0;
            if (if_trucncated == false)
                NUM_verts_j = 4;
            else
                NUM_verts_j = verts[j].NumVertsTruncated;

            float3 *Frac_verts_j = (float3 *)malloc(sizeof(float3) * NUM_verts_j); // = verts[j].verts_3D;

            if (if_trucncated == false)
                memcpy(Frac_verts_j, verts[j].Verts3D, sizeof(float3) * NUM_verts_j);
            else
                memcpy(Frac_verts_j, verts[j].Verts3DTruncated, sizeof(float3) * NUM_verts_j);

            for (int k = 0; k < NUM_verts_j; ++k)
            {
                Frac_verts_j[k] = make_float3(Frac_verts_j[k].x - Center_.x,
                                              Frac_verts_j[k].y - Center_.y,
                                              Frac_verts_j[k].z - Center_.z);
                Frac_verts_j[k] = cuDFNsys::ProductSquare3Float3(Mat_3_to_2,
                                                                 Frac_verts_j[k]);
            }

            //---------------start to identify
            float3 Intersection_xyPlane_Poly2[2];
            bool ik = cuDFNsys::Intersection3DPolyXYPlane(Frac_verts_j,
                                                          NUM_verts_j,
                                                          Intersection_xyPlane_Poly2,
                                                          _TOL_Intersection3DPolyXYPlane);
            if (ik == false)
            {
                free(Frac_verts);
                free(Frac_verts_j);
                free(Frac_verts_2D);
                continue;
            }
            float2 LineInf[2];
            LineInf[0] = make_float2(Intersection_xyPlane_Poly2[0].x, Intersection_xyPlane_Poly2[0].y);
            LineInf[1] = make_float2(Intersection_xyPlane_Poly2[1].x, Intersection_xyPlane_Poly2[1].y);

            float2 Intersection_Poly1_InfLine[2];
            bool jt = cuDFNsys::Intersection2DLine2DPoly(Frac_verts_2D,
                                                         NUM_verts,
                                                         LineInf,
                                                         Intersection_Poly1_InfLine,
                                                         _TOL_Intersection2DLine2DPoly);
            if (jt == false)
            {
                free(Frac_verts);
                free(Frac_verts_j);
                free(Frac_verts_2D);
                continue;
            }
            int sign_ = -1;
            float2 Intersection_two_segs[2];
            bool kh = cuDFNsys::IntersectionTwoCollinearSegs(LineInf,
                                                             Intersection_Poly1_InfLine,
                                                             Intersection_two_segs,
                                                             &sign_,
                                                             _TOL_IntersectionTwoCollinearSegs);

            if (kh == true)
            {
                if (sign_ == 2)
                {
                    float3 Intersection_f[2];
                    Intersection_f[0] = make_float3(Intersection_two_segs[0].x,
                                                    Intersection_two_segs[0].y, 0);
                    Intersection_f[1] = make_float3(Intersection_two_segs[1].x,
                                                    Intersection_two_segs[1].y, 0);

                    for (int k = 0; k < 2; ++k)
                    {
                        Intersection_f[k] = cuDFNsys::ProductSquare3Float3(Mat_2_to_3,
                                                                           Intersection_f[k]);
                        Intersection_f[k] = make_float3(Intersection_f[k].x + Center_.x, Intersection_f[k].y + Center_.y, Intersection_f[k].z + Center_.z);
                    }

                    pair<size_t, size_t> key_ = std::make_pair((size_t)i, (size_t)j);
                    pair<float3, float3> p;
                    p.first = Intersection_f[0];
                    p.second = Intersection_f[1];
                    pair<pair<size_t, size_t>, pair<float3, float3>> element_ =
                        std::make_pair(key_, p);
                    Intersection_map.insert(element_);

                    free(Frac_verts);
                    free(Frac_verts_j);
                    free(Frac_verts_2D);
                    continue;
                }
                else
                {

                    free(Frac_verts);
                    free(Frac_verts_j);
                    free(Frac_verts_2D);
                    continue;
                }
            }
            else
            {

                free(Frac_verts);
                free(Frac_verts_j);
                free(Frac_verts_2D);
                continue;
            }
        }
    }
}; // IdentifyIntersection
