#include "Fractures/TruncateFracture.cuh"

// ====================================================
// NAME:        TruncateFracture
// DESCRIPTION: Truncate a fracture in a DFN
// AUTHOR:      Tingchang YIN
// DATE:        06/04/2022
// ====================================================

__device__ __host__ bool cuDFNsys::TruncateFracture(cuDFNsys::Fracture *verts,
                                                    float L,
                                                    int plane,
                                                    int dir)
{
    bool if_touch = false;

    float bound = L * 0.5 * dir;
    //printf("bound %f\n", bound);

    float3 TT[8];

    int tmpcc = 0;

    for (int i = 0; i < verts->NumVertsTruncated; ++i)
    {
        float3 source = verts->Verts3DTruncated[i];
        float3 target = verts->Verts3DTruncated[(i + 1) % verts->NumVertsTruncated];

        //printf("The %d edge is being clipped:\n", i + 1);
        //printf("\tsource: %f, %f, %f\n", source.x, source.y, source.z);
        //printf("\ttarget: %f, %f, %f\n", target.x, target.y, target.z);

        float x1 = 0;
        float x2 = 0;
        if (plane == 0)
        {
            x1 = source.x;
            x2 = target.x;
        }
        else if (plane == 1)
        {
            x1 = source.y;
            x2 = target.y;
        }
        else if (plane == 2)
        {
            x1 = source.z;
            x2 = target.z;
        }

        float3 KK[2];
        int tmpll = -1;

        bool if_both_inside = false;
        bool if_both_outside = false;

        if (dir == -1) // sign = -1
        {
            if (x1 >= bound && x2 >= bound)
            {
                if_both_inside = true;
                if (x1 == bound || x2 == bound)
                    if_touch = true;
            }

            if (x1 < bound && x2 < bound)
                if_both_outside = true;
        }
        else if (dir == 1)
        {
            if (x1 <= bound && x2 <= bound)
            {
                if_both_inside = true;
                if (x1 == bound || x2 == bound)
                    if_touch = true;
            }

            if (x1 > bound && x2 > bound)
                if_both_outside = true;
        }

        //printf("both inside and both outside? %d %d\n", if_both_inside, if_both_outside);

        if (if_both_inside == true)
        {
            KK[0] = source;
            KK[1] = target;
            tmpll = 2;
        }
        else if (if_both_outside == true)
        {
            tmpll = 0;
        }
        else
        {
            float3 discardPNT, keepPNT;

            //------------------
            // if (abs(x1) > abs(bound))
            // {
            //     discardPNT = source;
            //     keepPNT = target;
            // }
            // else if (abs(x2) > abs(bound))
            // {
            //     discardPNT = target;
            //     keepPNT = source;
            // }
            //----------------------

            if (dir == -1) // sign = -1
            {
                if (x1 < bound && x2 >= bound)
                {
                    discardPNT = source;
                    keepPNT = target;
                }
                else if (x1 >= bound && x2 < bound)
                {
                    keepPNT = source;
                    discardPNT = target;
                }
            }
            else if (dir == 1)
            {
                if (x1 > bound && x2 <= bound)
                {
                    discardPNT = source;
                    keepPNT = target;
                }
                else if (x1 <= bound && x2 > bound)
                {
                    keepPNT = source;
                    discardPNT = target;
                }
            }

            //printf("Discard:\n\t%f, %f, %f\n", discardPNT.x, discardPNT.y, discardPNT.z);
            //printf("Keep:\n\t%f, %f, %f\n", keepPNT.x, keepPNT.y, keepPNT.z);
            // find intersection point
            // parametric function of a line
            float3 v = make_float3(discardPNT.x - keepPNT.x,
                                   discardPNT.y - keepPNT.y,
                                   discardPNT.z - keepPNT.z);

            float t = 0;
            float3 Intersection_ = make_float3(0, 0, 0);
            if (plane == 0)
            {
                t = (bound - keepPNT.x) / v.x;
                Intersection_.x = bound;
                Intersection_.y = keepPNT.y + t * v.y;
                Intersection_.z = keepPNT.z + t * v.z;
            }
            else if (plane == 1)
            {
                t = (bound - keepPNT.y) / v.y;
                Intersection_.x = keepPNT.x + t * v.x;
                Intersection_.y = bound;
                Intersection_.z = keepPNT.z + t * v.z;
            }
            else if (plane == 2)
            {
                t = (bound - keepPNT.z) / v.z;
                Intersection_.x = keepPNT.x + t * v.x;
                Intersection_.y = keepPNT.y + t * v.y;
                Intersection_.z = bound;
            }

            /*if (abs(x1) > abs(bound))
            {
                KK[0] = Intersection_;
                KK[1] = keepPNT;
            }
            else if (abs(x2) > abs(bound))
            {
                KK[1] = Intersection_;
                KK[0] = keepPNT;
            }*/

            if (dir == -1)
            {
                if (x1 < bound)
                {
                    KK[0] = Intersection_;
                    KK[1] = keepPNT;
                }
                else if (x2 < bound)
                {
                    KK[1] = Intersection_;
                    KK[0] = keepPNT;
                }
            }
            else if (dir == 1)
            {
                if (x1 > bound)
                {
                    KK[0] = Intersection_;
                    KK[1] = keepPNT;
                }
                else if (x2 > bound)
                {
                    KK[1] = Intersection_;
                    KK[0] = keepPNT;
                }
            }

            tmpll = 2;
            if_touch = true;
        }

        if (tmpll != 0)
        {
            //printf("inserting pnts\n");
            for (int j = 0; j < tmpll; ++j)
            {
                float3 PNT_this = KK[j];

                bool if_duplicated = false;

                for (int k = 0; k < tmpcc; ++k)
                {
                    float3 PNT_e = TT[k];

                    float3 dd = make_float3(PNT_this.x - PNT_e.x, PNT_this.y - PNT_e.y, PNT_this.z - PNT_e.z);
                    float distance = pow(dd.x * dd.x + dd.y * dd.y + dd.z * dd.z, 0.5);
                    if (distance == 0)
                    {
                        if_duplicated = true;
                        break;
                    }
                }

                if (if_duplicated == false)
                {
                    TT[tmpcc] = PNT_this;
                    //printf("\t%f, %f, %f\n", PNT_this.x, PNT_this.y, PNT_this.z);
                    tmpcc++;
                }
            }
            //printf("------\n\n");
        }
    }

    verts->NumVertsTruncated = tmpcc;

    //printf("[");
    for (int i = 0; i < tmpcc; ++i)
    {
        verts->Verts3DTruncated[i] = TT[i];
        //printf("%f, %f, %f;\n", TT[i].x, TT[i].y, TT[i].z);
    };
    //printf("]");
    return if_touch;
};