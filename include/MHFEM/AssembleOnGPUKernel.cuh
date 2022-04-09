///////////////////////////////////////////////////////////////////
// NAME:              AssembleOnGPUKernel.cuh
//
// PURPOSE:           kernel function to assemble matrix
//
// FUNCTIONS/OBJECTS: AssembleOnGPUKernel
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../GlobalDef/GlobalDef.cuh"
#include "../Mesh/Mesh.cuh"
#include "Triplet.cuh"
#include "StimaA.cuh"

namespace cuDFNsys
{
__global__ void AssembleOnGPUKernel(cuDFNsys::Triplet *tri_dev,
                                    cuDFNsys::EleCoor *coord_2D_dev,
                                    cuDFNsys::EleEdgeAttri *Edge_attri,
                                    float *Conduc_Frac_dev,
                                    //int *neuman_sep_dev,
                                    int NUM_sep_edges,
                                    int NUM_eles,
                                    int NUM_glob_interior_edges,
                                    float P_in,
                                    float P_out);
}; // namespace cuDFNsys