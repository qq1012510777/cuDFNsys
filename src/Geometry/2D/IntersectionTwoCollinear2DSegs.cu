#include "Geometry/2D/IntersectionTwoCollinear2DSegs.cuh"

// ====================================================
// NAME:        IntersectionTwoCollinear2DSegs
// DESCRIPTION: Identify intersection between two collinear segments
// AUTHOR:      Tingchang YIN
// DATE:        07/04/2022
// ====================================================
template <typename T>
__device__ __host__ bool cuDFNsys::IntersectionTwoCollinear2DSegs(cuDFNsys::Vector2<T> *Seg_1,
                                                                cuDFNsys::Vector2<T> *Seg_2,
                                                                cuDFNsys::Vector2<T> *intersection,
                                                                int *sign,
                                                                T _TOL_)
{
    cuDFNsys::Vector2<T> A, B, C, D;

    cuDFNsys::Vector2<T> v1 = cuDFNsys::MakeVector2<T>(Seg_1[0].x - Seg_1[1].x,
                                                       Seg_1[0].y - Seg_1[1].y);
    cuDFNsys::Vector2<T> v2 = cuDFNsys::MakeVector2<T>(Seg_2[0].x - Seg_2[1].x,
                                                       Seg_2[0].y - Seg_2[1].y);

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
        T pb_pc = B.x - C.x;
        T pd_pa = D.x - A.x;

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

            cuDFNsys::Vector2<T> dd = cuDFNsys::MakeVector2<T>(intersection[0].x - intersection[1].x, intersection[0].y - intersection[1].y);
            T distance = pow(dd.x * dd.x + dd.y * dd.y, 0.5);

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
        T pb_pc = B.y - C.y;
        T pd_pa = D.y - A.y;

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

            cuDFNsys::Vector2<T> dd = cuDFNsys::MakeVector2<T>(intersection[0].x - intersection[1].x, intersection[0].y - intersection[1].y);
            T distance = pow(dd.x * dd.x + dd.y * dd.y, 0.5);

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
}; // IntersectionTwoCollinear2DSegs
template __device__ __host__ bool cuDFNsys::IntersectionTwoCollinear2DSegs<double>(cuDFNsys::Vector2<double> *Seg_1,
                                                                                 cuDFNsys::Vector2<double> *Seg_2,
                                                                                 cuDFNsys::Vector2<double> *intersection,
                                                                                 int *sign,
                                                                                 double _TOL_);
template __device__ __host__ bool cuDFNsys::IntersectionTwoCollinear2DSegs<float>(cuDFNsys::Vector2<float> *Seg_1,
                                                                                cuDFNsys::Vector2<float> *Seg_2,
                                                                                cuDFNsys::Vector2<float> *intersection,
                                                                                int *sign,
                                                                                float _TOL_);