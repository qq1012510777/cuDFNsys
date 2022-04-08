#include "Geometry/2D/Intersection2DLine2DPoly.cuh"

// ====================================================
// NAME:        Intersection2DLine2DPoly
// DESCRIPTION: Identify Intersections between
//              2D Line (infinite) and 2D Polygon
// AUTHOR:      Tingchang YIN
// DATE:        07/04/2022
// ====================================================
__device__ __host__ bool cuDFNsys::Intersection2DLine2DPoly(float2 *Poly2D,
                                                            int NUM_verts,
                                                            float2 *Line,
                                                            float2 *intersection_k,
                                                            float _TOL_)
{
    int tmpcc = 0;
    float2 tmp_pnts[10];
    for (int i = 0; i < NUM_verts; ++i)
    {
        float2 Seg[2], intersection[2];
        Seg[0] = Poly2D[i];
        Seg[1] = Poly2D[(i + 1) % NUM_verts];

        int sign__ = 0;

        bool ik = cuDFNsys::Intersection2DLine2DSeg(Line,
                                                    Seg,
                                                    &sign__,
                                                    intersection,
                                                    _TOL_);

        if (ik == true)
        {

            if (sign__ == 1)
            {

                bool if_duplicated = false;
                for (int j = 0; j < tmpcc; ++j)
                {
                    float2 dd = make_float2(intersection[0].x - tmp_pnts[j].x, intersection[0].y - tmp_pnts[j].y);
                    float distance = pow(dd.x * dd.x + dd.y * dd.y, 0.5);

                    if (distance < _TOL_)
                    {
                        if_duplicated = true;
                        break;
                    }
                }

                if (if_duplicated == false)
                {
                    tmp_pnts[tmpcc] = intersection[0];
                    tmpcc++;
                }
            }
        }
    }

    if (tmpcc == 2)
    {
        intersection_k[0] = tmp_pnts[0];
        intersection_k[1] = tmp_pnts[1];
        return true;
    }

    if (tmpcc == 1)
    {
        return false;
    }

    if (tmpcc == 0)
    {
        return false;
    }

    if (tmpcc > 2)
    {
        //printf("Error! An infinite line (2D) cannot intersect a 2D polygon with more than two points!\n");
        //exit(0);

        uint2 pair__;

        float dist = 0;
        for (int k = 0; k < tmpcc - 1; ++k)
        {
            for (int h = k + 1; h < tmpcc; ++h)
            {
                float2 kh = make_float2(tmp_pnts[k].x - tmp_pnts[h].x, tmp_pnts[k].y - tmp_pnts[h].y);
                float distrr = sqrt(kh.x * kh.x + kh.y * kh.y);

                if (distrr > dist)
                {
                    dist = distrr;
                    pair__.x = k;
                    pair__.y = h;
                }
            }
        };

        intersection_k[0] = tmp_pnts[pair__.x];
        intersection_k[1] = tmp_pnts[pair__.y];
        return true;
    }

    return false;
}; // Intersection2DLine2DPoly