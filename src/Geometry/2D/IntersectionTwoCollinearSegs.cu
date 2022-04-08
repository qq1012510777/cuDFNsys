#include "Geometry/2D/IntersectionTwoCollinearSegs.cuh"

// ====================================================
// NAME:        IntersectionTwoCollinearSegs
// DESCRIPTION: Identify intersection between two collinear segments
// AUTHOR:      Tingchang YIN
// DATE:        07/04/2022
// ====================================================

__device__ __host__ bool cuDFNsys::IntersectionTwoCollinearSegs(float2 *Seg_1,
                                                                float2 *Seg_2,
                                                                float2 *intersection,
                                                                int *sign,
                                                                float _TOL_)
{
    float2 A, B, C, D;

    float2 v1 = make_float2(Seg_1[0].x - Seg_1[1].x, Seg_1[0].y - Seg_1[1].y);
    float2 v2 = make_float2(Seg_2[0].x - Seg_2[1].x, Seg_2[0].y - Seg_2[1].y);

    if (abs(v1.x) > _TOL_)
    {
        if (v1.x > 0)
        {
            A = Seg_1[1];
            B = Seg_1[0];
        }
        else
        {
            A = Seg_1[0];
            B = Seg_1[1];
        }

        if (v2.x > 0)
        {
            C = Seg_2[1];
            D = Seg_2[0];
        }
        else
        {
            C = Seg_2[0];
            D = Seg_2[1];
        }
        //--------------------
        float pb_pc = B.x - C.x;
        float pd_pa = D.x - A.x;

        if (pb_pc >= 0 && pd_pa >= 0)
        {
            if (A.x > C.x)
                intersection[0] = A;
            else
                intersection[0] = C;

            if (B.x < D.x)
                intersection[1] = B;
            else
                intersection[1] = D;

            float2 dd = make_float2(intersection[0].x - intersection[1].x, intersection[0].y - intersection[1].y);
            float distance = pow(dd.x * dd.x + dd.y * dd.y, 0.5);

            if (distance < _TOL_)
            {
                *sign = 1;
            }
            else
                *sign = 2;

            return true;
        }
        else
        {
            *sign = -1;
            return false;
        }
    }
    else if (abs(v1.y) > _TOL_)
    {
        if (v1.y > 0)
        {
            A = Seg_1[1];
            B = Seg_1[0];
        }
        else
        {
            A = Seg_1[0];
            B = Seg_1[1];
        }

        if (v2.y > 0)
        {
            C = Seg_2[1];
            D = Seg_2[0];
        }
        else
        {
            C = Seg_2[0];
            D = Seg_2[1];
        }

        //--------------------
        float pb_pc = B.y - C.y;
        float pd_pa = D.y - A.y;

        if (pb_pc >= 0 && pd_pa >= 0)
        {
            if (A.y > C.y)
                intersection[0] = A;
            else
                intersection[0] = C;

            if (B.y < D.y)
                intersection[1] = B;
            else
                intersection[1] = D;

            float2 dd = make_float2(intersection[0].x - intersection[1].x, intersection[0].y - intersection[1].y);
            float distance = pow(dd.x * dd.x + dd.y * dd.y, 0.5);

            if (distance < _TOL_)
            {
                *sign = 1;
            }
            else
                *sign = 2;
            return true;
        }
        else
        {
            *sign = -1;
            return false;
        }
    }
    else
    {
        *sign = -1;
        // printf("Error! Two segments are two short!\n");
        return false;
    }

    *sign = -1;
    return false;
}; // IntersectionTwoCollinearSegs