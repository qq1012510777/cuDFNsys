#include "MHFEM/AssembleOnGPUKernel.cuh"

// ====================================================
// NAME:        AssembleOnGPUKernel
// DESCRIPTION: global function to assemble matrix
// AUTHOR:      Tingchang YIN
// DATE:        09/04/2022
// ====================================================
template <typename T>
__global__ void cuDFNsys::AssembleOnGPUKernel(cuDFNsys::Triplet<T> *tri_dev,
                                              cuDFNsys::EleCoor<T> *coord_2D_dev,
                                              cuDFNsys::EleEdgeAttri *Edge_attri,
                                              T *Conduc_Frac_dev,
                                              //int *neuman_sep_dev,
                                              int NUM_sep_edges,
                                              int NUM_eles,
                                              int NUM_glob_interior_edges,
                                              T P_in,
                                              T P_out)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= NUM_eles)
        return;
    int I[3] = {i * 3 + 1, // 2
                i * 3 + 2, // 3
                i * 3};    // 1

    cuDFNsys::EleCoor<T> coord = coord_2D_dev[i];

    T A[3][3];
    cuDFNsys::StimaA<T>(coord, A);

    size_t j = i * 21;
    size_t j_tt = j;

    for (size_t ik = 0; ik < 3; ++ik)
        for (size_t jk = 0; jk < 3; ++jk)
        {
            tri_dev[j].row = I[ik];
            tri_dev[j].col = I[jk];
            tri_dev[j].val = -1.0f / Conduc_Frac_dev[i] * A[ik][jk];

            int edge_1 = tri_dev[j].row % 3;
            int edge_2 = tri_dev[j].col % 3;
            if (Edge_attri[i].e[edge_1] == 2 && Edge_attri[i].e[edge_2] == 2 && edge_1 == edge_2)
                tri_dev[j].val = 1.0;
            else if (edge_1 != edge_2 && (Edge_attri[i].e[edge_1] == 2 || Edge_attri[i].e[edge_2] == 2))
                tri_dev[j].row = -1;

            j++;
        }

    T B[3] = {0};
    B[0] = pow(pow(coord.x[2] - coord.x[1], 2) + pow(coord.y[2] - coord.y[1], 2), 0.5);
    B[1] = pow(pow(coord.x[0] - coord.x[2], 2) + pow(coord.y[0] - coord.y[2], 2), 0.5);
    B[2] = pow(pow(coord.x[1] - coord.x[0], 2) + pow(coord.y[1] - coord.y[0], 2), 0.5);

    for (size_t ik = 0; ik < 3; ++ik)
    {
        tri_dev[j].row = I[ik];
        tri_dev[j].col = i + NUM_sep_edges;
        tri_dev[j].val = B[ik];
        j++;
        tri_dev[j].row = i + NUM_sep_edges;
        tri_dev[j].col = I[ik];
        tri_dev[j].val = B[ik];
        j++;

        int edge_1 = tri_dev[j - 2].row % 3;
        if (Edge_attri[i].e[edge_1] == 2)
        {
            tri_dev[j - 2].row = -1;
            tri_dev[j - 1].row = -1;
        }
    }

    T P_in_out[2] = {P_in, P_out};

    for (size_t ik = 0; ik < 3; ++ik)
    {
        int ek = Edge_attri[i].e[ik];
        int NO_ = Edge_attri[i].no[ik] - 1;

        if (ek == 3)
        {
            tri_dev[j].row = NO_ + NUM_sep_edges + NUM_eles;
            tri_dev[j].col = I[(ik + 2) % 3];
            tri_dev[j].val = -B[(ik + 2) % 3];
            j++;
            tri_dev[j].col = NO_ + NUM_sep_edges + NUM_eles;
            tri_dev[j].row = I[(ik + 2) % 3];
            tri_dev[j].val = -B[(ik + 2) % 3];
            j++;
        }
        else if (ek == 0 || ek == 1) // in or out
        {
            tri_dev[j].row = NO_;
            tri_dev[j].col = NUM_sep_edges +
                             NUM_eles +
                             NUM_glob_interior_edges + 2;
            tri_dev[j].val = P_in_out[ek] * B[(ik + 2) % 3];
            j++;
        }
        //else if (ek == 2) // neumann
        //neuman_sep_dev[NO_] = 1 + NO_;
    };

    for (size_t k = j; k < j_tt + 21; ++k)
        tri_dev[k].row = -1;
}; // AssembleOnGPUKernel
template __global__ void cuDFNsys::AssembleOnGPUKernel<double>(cuDFNsys::Triplet<double> *tri_dev,
                                                               cuDFNsys::EleCoor<double> *coord_2D_dev,
                                                               cuDFNsys::EleEdgeAttri *Edge_attri,
                                                               double *Conduc_Frac_dev,
                                                               //int *neuman_sep_dev,
                                                               int NUM_sep_edges,
                                                               int NUM_eles,
                                                               int NUM_glob_interior_edges,
                                                               double P_in,
                                                               double P_out);
template __global__ void cuDFNsys::AssembleOnGPUKernel<float>(cuDFNsys::Triplet<float> *tri_dev,
                                                              cuDFNsys::EleCoor<float> *coord_2D_dev,
                                                              cuDFNsys::EleEdgeAttri *Edge_attri,
                                                              float *Conduc_Frac_dev,
                                                              //int *neuman_sep_dev,
                                                              int NUM_sep_edges,
                                                              int NUM_eles,
                                                              int NUM_glob_interior_edges,
                                                              float P_in,
                                                              float P_out);