#include "MHFEM/StimaA.cuh"

// ====================================================
// NAME:        StimaA
// DESCRIPTION: Assemble submatrix A for each element
// AUTHOR:      Tingchang YIN
// DATE:        09/04/2022
// ====================================================
__device__ __host__ void cuDFNsys::StimaA(cuDFNsys::EleCoor coord, float A[3][3])
{
    float P1_P2[2], P2_P1[2], P3_P2[2], P2_P3[2], P3_P1[2], P1_P3[2];

    P1_P2[0] = coord.x[0] - coord.x[1];
    P1_P2[1] = coord.y[0] - coord.y[1];

    P3_P2[0] = coord.x[2] - coord.x[1];
    P3_P2[1] = coord.y[2] - coord.y[1];

    P3_P1[0] = coord.x[2] - coord.x[0];
    P3_P1[1] = coord.y[2] - coord.y[0];

    P2_P1[0] = -P1_P2[0];
    P2_P1[1] = -P1_P2[1];

    P2_P3[0] = -P3_P2[0];
    P2_P3[1] = -P3_P2[1];

    P1_P3[0] = -P3_P1[0];
    P1_P3[1] = -P3_P1[1];

    float N_trs[3][6] = {0, 0, P2_P1[0], P2_P1[1], P3_P1[0], P3_P1[1],
                         P1_P2[0], P1_P2[1], 0, 0, P3_P2[0], P3_P2[1],
                         P1_P3[0], P1_P3[1], P2_P3[0], P2_P3[1], 0, 0};

    float N[6][3] = {0, P1_P2[0], P1_P3[0],
                     0, P1_P2[1], P1_P3[1],
                     P2_P1[0], 0, P2_P3[0],
                     P2_P1[1], 0, P2_P3[1],
                     P3_P1[0], P3_P2[0], 0,
                     P3_P1[1], P3_P2[1], 0};

    float M[6][6] = {
        2, 0, 1, 0, 1, 0,
        0, 2, 0, 1, 0, 1,
        1, 0, 2, 0, 1, 0,
        0, 1, 0, 2, 0, 1,
        1, 0, 1, 0, 2, 0,
        0, 1, 0, 1, 0, 2};

    float T[3][3] = {coord.x[0], coord.x[1], coord.x[2],
                     coord.y[0], coord.y[1], coord.y[2],
                     1, 1, 1};

    float D[3][3] = {sqrt(P3_P2[0] * P3_P2[0] + P3_P2[1] * P3_P2[1]), 0, 0,
                     0, sqrt(P1_P3[0] * P1_P3[0] + P1_P3[1] * P1_P3[1]), 0,
                     0, 0, sqrt(P1_P2[0] * P1_P2[0] + P1_P2[1] * P1_P2[1])};

    //-----
    float Mat_1[3][6] = {0};
    for (size_t ik = 0; ik < 3; ik++)
        for (size_t jk = 0; jk < 6; jk++)
        {
            for (size_t k = 0; k < 3; k++)
                Mat_1[ik][jk] += D[ik][k] * N_trs[k][jk];
        }
    //-------------------
    float Mat_2[3][6] = {0};
    for (size_t ik = 0; ik < 3; ik++)
        for (size_t jk = 0; jk < 6; jk++)
        {
            for (size_t k = 0; k < 6; k++)
                Mat_2[ik][jk] += Mat_1[ik][k] * M[k][jk];
        }
    //----------

    float Mat_3[3][3] = {0};
    //Mat_multiply(Mat_2, 3, 6, N, 3, Mat_3);
    for (size_t ik = 0; ik < 3; ik++)
        for (size_t jk = 0; jk < 3; jk++)
        {
            for (size_t k = 0; k < 6; k++)
                Mat_3[ik][jk] += Mat_2[ik][k] * N[k][jk];
        }
    //-------

    float Mat_4[3][3] = {0};
    //Mat_multiply(Mat_3, 3, 3, D, 3, Mat_4);
    for (size_t ik = 0; ik < 3; ik++)
        for (size_t jk = 0; jk < 3; jk++)
        {
            for (size_t k = 0; k < 3; k++)
                Mat_4[ik][jk] += Mat_3[ik][k] * D[k][jk];
        }
    //-------

    float determinant_ = T[0][0] * ((T[1][1] * T[2][2]) - (T[2][1] * T[1][2])) - T[0][1] * (T[1][0] * T[2][2] - T[2][0] * T[1][2]) + T[0][2] * (T[1][0] * T[2][1] - T[2][0] * T[1][1]);

    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
            A[i][j] = Mat_4[i][j] / (24 * determinant_);
}; // StimaA