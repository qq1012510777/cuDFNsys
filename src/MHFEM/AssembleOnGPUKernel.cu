/****************************************************************************
* cuDFNsys - simulating flow and transport in 3D fracture networks          *
* Copyright (C) 2022, Tingchang YIN, Sergio GALINDO-TORRES                  *
*                                                                           *
* This program is free software: you can redistribute it and/or modify      *
* it under the terms of the GNU Affero General Public License as            *
* published by the Free Software Foundation, either version 3 of the        *
* License, or (at your option) any later version.                           *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU Affero General Public License for more details.                       *
*                                                                           *
* You should have received a copy of the GNU Affero General Public License  *
* along with this program.  If not, see <https://www.gnu.org/licenses/>.    *
*****************************************************************************/

#include "MHFEM/AssembleOnGPUKernel.cuh"

// ====================================================
// NAME:        AssembleOnGPUKernel
// DESCRIPTION: global function to assemble matrix
// AUTHOR:      Tingchang YIN
// DATE:        09/04/2022
// ====================================================
template <typename T>
__global__ void cuDFNsys::AssembleOnGPUKernel(
    cuDFNsys::Triplet<T> *tri_dev, cuDFNsys::EleCoor<T> *coord_2D_dev,
    cuDFNsys::EleEdgeAttri *Edge_attri, T *Conduc_Frac_dev,
    //int *neuman_sep_dev,
    int NUM_sep_edges, int NUM_eles, int Dim, T P_in, T P_out, bool IfPeriodic,
    T mu_over_RhoGravity, T const_q)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= NUM_eles)
        return;
    int I[3] = {i * 3 + 1,  // 2
                i * 3 + 2,  // 3
                i * 3 + 0}; // 1

    cuDFNsys::EleCoor<T> coord = coord_2D_dev[i];

    T A[3][3];
    cuDFNsys::StimaA<T>(coord, A);

    // if (i == 0)
    // {
    //     printf("%f, %f, %f\n%f, %f, %f\n%f, %f, %f\n***\n", A[0][0], A[0][1],
    //            A[0][2], A[1][0], A[1][1], A[1][2], A[2][0], A[2][1], A[2][2]);
    // }

    size_t j = i * 18;
    size_t j_tt = j;

    for (size_t ik = 0; ik < 3; ++ik)
        for (size_t jk = 0; jk < 3; ++jk)
        {
            tri_dev[j].row = I[ik];
            tri_dev[j].col = I[jk];
            tri_dev[j].val =
                (T)-1.0 / (Conduc_Frac_dev[i] * mu_over_RhoGravity) * A[ik][jk];

            // int edge_1 = tri_dev[j].row % 3;
            // int edge_2 = tri_dev[j].col % 3;
            // if (Edge_attri[i].e[edge_1] == 2 && Edge_attri[i].e[edge_2] == 2 && edge_1 == edge_2) // non-flux ???
            //     tri_dev[j].val = 1.0;
            // else if (edge_1 != edge_2 && (Edge_attri[i].e[edge_1] == 2 || Edge_attri[i].e[edge_2] == 2))
            //     tri_dev[j].row = -1;

            j++;
        }

    T B[3] = {0};
    B[0] = pow(
        pow(coord.x[2] - coord.x[1], 2) + pow(coord.y[2] - coord.y[1], 2), 0.5);
    B[1] = pow(
        pow(coord.x[0] - coord.x[2], 2) + pow(coord.y[0] - coord.y[2], 2), 0.5);
    B[2] = pow(
        pow(coord.x[1] - coord.x[0], 2) + pow(coord.y[1] - coord.y[0], 2), 0.5);

    for (size_t ik = 0; ik < 3; ++ik)
    {
        /// tri_dev[j].row = I[ik];
        /// tri_dev[j].col = i + NUM_sep_edges; // element NO
        /// tri_dev[j].val = B[ik];
        /// j++;
        tri_dev[j].row = i + NUM_sep_edges;
        tri_dev[j].col = I[ik];
        tri_dev[j].val = B[ik];
        j++;

        /// int edge_1 = tri_dev[j - 2].row % 3;
        /// if (Edge_attri[i].e[edge_1] == 2)
        /// {
        ///     tri_dev[j - 2].row = -1;
        ///     tri_dev[j - 1].row = -1;
        /// }
    }

    T P_in_out[2] = {P_in, P_out};

    for (size_t ik = 0; ik < 3; ++ik)
    {
        int ek = Edge_attri[i].e[I[ik] % 3];
        int NO_ =
            Edge_attri[i]
                .no[I[ik] % 3]; // global edge NO for all edges, no separating
        // printf("NO: %d, element %d\n", NO_, i);

        if (ek == 3)
        {
            tri_dev[j].row = NO_ + NUM_sep_edges + NUM_eles;
            tri_dev[j].col = I[ik];
            tri_dev[j].val = -B[ik]; // -B[(ik + 2) % 3];
            j++;
            /// tri_dev[j].col = NO_ + NUM_sep_edges + NUM_eles;
            /// tri_dev[j].row = I[ik];  //I[(ik + 2) % 3];
            /// tri_dev[j].val = -B[ik]; //-B[(ik + 2) % 3];
            /// j++;
        }
        else if (ek == 0 || ek == 1) // in or out
        {
            if (!IfPeriodic)
            {
                tri_dev[j].row = I[ik]; // I[(ik + 2) % 3];
                tri_dev[j].col = Dim + 2;
                tri_dev[j].val = P_in_out[ek] * B[ik]; //B[(ik + 2) % 3];
                j++;
            }
            else
            {
                tri_dev[j].row = NO_ + NUM_sep_edges + NUM_eles;
                tri_dev[j].col = I[ik];  // I[(ik + 2) % 3];
                tri_dev[j].val = -B[ik]; //-B[(ik + 2) % 3];
                j++;

                tri_dev[j].row = NO_ + NUM_sep_edges + NUM_eles;
                tri_dev[j].col = Dim + 2;
                tri_dev[j].val =
                    (ek == 0 ? -1 : 1) * const_q * -B[ik]; //B[(ik + 2) % 3];
                j++;
            }
        }
        else if (ek == 2) // neumann
        {
            // neuman_sep_dev[NO_] = 1 + NO_;
            tri_dev[j].row = NO_ + NUM_sep_edges + NUM_eles;
            tri_dev[j].col = I[ik];  // I[(ik + 2) % 3];
            tri_dev[j].val = -B[ik]; //-B[(ik + 2) % 3];
            j++;

            /// tri_dev[j].col = NO_ + NUM_sep_edges + NUM_eles;
            /// tri_dev[j].row = I[ik];  // I[(ik + 2) % 3];
            /// tri_dev[j].val = -B[ik]; //-B[(ik + 2) % 3];
            /// j++;

            tri_dev[j].row = NO_ + NUM_sep_edges + NUM_eles;
            tri_dev[j].col = Dim + 2;
            tri_dev[j].val = 0 * -B[ik]; //B[(ik + 2) % 3];
            j++;
        }
    };

    for (size_t k = j; k < j_tt + 18; ++k)
        tri_dev[k].row = -1;
}; // AssembleOnGPUKernel
template __global__ void cuDFNsys::AssembleOnGPUKernel<double>(
    cuDFNsys::Triplet<double> *tri_dev, cuDFNsys::EleCoor<double> *coord_2D_dev,
    cuDFNsys::EleEdgeAttri *Edge_attri, double *Conduc_Frac_dev,
    //int *neuman_sep_dev,
    int NUM_sep_edges, int NUM_eles, int Dim, double P_in, double P_out,
    bool IfPeriodic, double mu_over_RhoGravity, double const_q);
template __global__ void cuDFNsys::AssembleOnGPUKernel<float>(
    cuDFNsys::Triplet<float> *tri_dev, cuDFNsys::EleCoor<float> *coord_2D_dev,
    cuDFNsys::EleEdgeAttri *Edge_attri, float *Conduc_Frac_dev,
    //int *neuman_sep_dev,
    int NUM_sep_edges, int NUM_eles, int Dim, float P_in, float P_out,
    bool IfPeriodic, float mu_over_RhoGravity, float const_q);