#include "Geometry/2D/Intersection2DLine2DSeg.cuh"

// ====================================================
// NAME:        Intersection2DLine2DSeg
// DESCRIPTION: Identify Intersections between
//              2D line (infinite) and 2D segment
// AUTHOR:      Tingchang YIN
// DATE:        07/04/2022
// ====================================================

__device__ __host__ bool cuDFNsys::Intersection2DLine2DSeg(float2 *Line,
                                                           float2 *Seg,
                                                           int *sign_, // 1, pnt; 2, seg; 3, none
                                                           float2 *intersection,
                                                           float _TOL_)
{
    float2 directional_v = make_float2(Line[0].x - Line[1].x, Line[0].y - Line[1].y);

    float2 Normal_To_Line = make_float2(-1.0 * directional_v.y, directional_v.x);

    float2 p1 = Seg[0];
    float2 p2 = Seg[1];

    float2 p = make_float2((Line[0].x + Line[1].x) * 0.5, (Line[0].y + Line[1].y) * 0.5);

    float2 j1 = make_float2(p1.x - p.x, p1.y - p.y);
    float2 j2 = make_float2(p2.x - p.x, p2.y - p.y);

    float g1 = Normal_To_Line.x * j1.x + Normal_To_Line.y * j1.y;
    float g2 = Normal_To_Line.x * j2.x + Normal_To_Line.y * j2.y;

    if (abs(g1) < _TOL_ && abs(g2) < _TOL_) // L is containing the Seg
    {
        *sign_ = 2;
        intersection[0] = Seg[0];
        intersection[1] = Seg[1];
        return true;
    }
    else if (abs(g1) < _TOL_ && abs(g2) > _TOL_) // one point lies on the L
    {
        *sign_ = 1;
        intersection[0] = Seg[0];
        return true;
    }
    else if (abs(g1) > _TOL_ && abs(g2) < _TOL_) // another point lies on the L
    {
        *sign_ = 1;
        intersection[0] = Seg[1];
        return true;
    }

    if ((g1 / abs(g1)) + (g2 / abs(g2)) == 0) // oppsite signs
    {
        *sign_ = 1;
        // now intersection between two infinite lines is the intersection we need
        float a1 = Line[0].x;
        float a2 = Line[0].y;
        float b1 = Line[1].x;
        float b2 = Line[1].y;

        float c1 = Seg[0].x;
        float c2 = Seg[0].y;
        float d1 = Seg[1].x;
        float d2 = Seg[1].y;

        float A1 = a2 - b2;
        float B1 = b1 - a1;
        float C1 = (A1 * a1) + (B1 * a2);
        float A2 = c2 - d2;
        float B2 = d1 - c1;
        float C2 = (A2 * c1) + (B2 * c2);

        float x = (C1 * B2 - C2 * B1) / (A1 * B2 - A2 * B1);
        float y = (C2 * A1 - C1 * A2) / (A1 * B2 - B1 * A2);
        intersection[0] = make_float2(x, y);
        return true;
    }

    *sign_ = -1;
    return false;
}; // Intersection2DLine2DSeg