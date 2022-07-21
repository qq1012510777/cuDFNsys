#include "Geometry/3D/Intersection3DPolyXYPlane.cuh"

// ====================================================
// NAME:        Intersection3DPolyXYPlane
// DESCRIPTION: Identify intersection between
//              a 3D polygon and the XY plane
// AUTHOR:      Tingchang YIN
// DATE:        07/04/2022
// ====================================================
template <typename T>
__device__ __host__ bool cuDFNsys::Intersection3DPolyXYPlane(cuDFNsys::Vector3<T> *Poly,
                                                             int NUM_vert,
                                                             cuDFNsys::Vector3<T> *Intersection,
                                                             T _TOL_)
{
    cuDFNsys::Vector3<T> pnt_t[10];
    int tmp_c = 0;

    for (int i = 0; i < NUM_vert; ++i)
    {
        cuDFNsys::Vector3<T> Seg[2], Intersec_PNT[2];
        Seg[0] = Poly[i];
        Seg[1] = Poly[(i + 1) % NUM_vert];

        int sign_;

        bool ik = cuDFNsys::Intersection3DSegXYPlane<T>(Seg, Intersec_PNT, &sign_, _TOL_);

        if (sign_ == 1)
        {
            bool if_duplicated = false;

            for (int j = 0; j < tmp_c; ++j)
            {
                cuDFNsys::Vector3<T> dd = cuDFNsys::MakeVector3(Intersec_PNT[0].x - pnt_t[j].x,
                                                                Intersec_PNT[0].y - pnt_t[j].y,
                                                                Intersec_PNT[0].z - pnt_t[j].z);

                T dis_value = pow(pow(dd.x, 2) + pow(dd.y, 2) + pow(dd.z, 2), 0.5);

                if (dis_value < _TOL_)
                {
                    if_duplicated = true;
                    continue;
                }
            }

            if (if_duplicated == false)
            {
                pnt_t[tmp_c] = Intersec_PNT[0];
                tmp_c++;
            }
        }
    }

    if (tmp_c == 2)
    {
        Intersection[0] = pnt_t[0];
        Intersection[1] = pnt_t[1];
        return true;
    }

    if (tmp_c == 0)
        return false;

    if (tmp_c == 1)
        return false;

    if (tmp_c > 2)
    {
        //printf("error, a 3D polygon cannot intersect the xy plane with more than two points!\n");
        //exit(0);
        uint2 pair__;

        T dist = 0;
        for (int k = 0; k < tmp_c - 1; ++k)
        {
            for (int h = k + 1; h < tmp_c; ++h)
            {
                cuDFNsys::Vector3<T> kh = cuDFNsys::MakeVector3(pnt_t[k].x - pnt_t[h].x, pnt_t[k].y - pnt_t[h].y, pnt_t[k].z - pnt_t[h].z);
                T distrr = sqrt(kh.x * kh.x + kh.y * kh.y + kh.z * kh.z);

                if (distrr > dist)
                {
                    dist = distrr;
                    pair__.x = k;
                    pair__.y = h;
                }
            }
        };

        Intersection[0] = pnt_t[pair__.x];
        Intersection[1] = pnt_t[pair__.y];
        return true;
    }

    return false;
};
template __device__ __host__ bool cuDFNsys::Intersection3DPolyXYPlane<double>(cuDFNsys::Vector3<double> *Poly,
                                                                              int NUM_vert,
                                                                              cuDFNsys::Vector3<double> *Intersection,
                                                                              double _TOL_);
template __device__ __host__ bool cuDFNsys::Intersection3DPolyXYPlane<float>(cuDFNsys::Vector3<float> *Poly,
                                                                             int NUM_vert,
                                                                             cuDFNsys::Vector3<float> *Intersection,
                                                                             float _TOL_);