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
#include "../DataTypeSelector/DataTypeSelector.cuh"
#include "../GlobalDef/GlobalDef.cuh"
#include "../Mesh/Mesh.cuh"
#include "StimaA.cuh"
#include "Triplet.cuh"

namespace cuDFNsys
{
template <typename T>
__global__ void AssembleOnGPUKernel(cuDFNsys::Triplet<T> *tri_dev,
                                    cuDFNsys::EleCoor<T> *coord_2D_dev,
                                    cuDFNsys::EleEdgeAttri *Edge_attri,
                                    T *Conduc_Frac_dev,
                                    //int *neuman_sep_dev,
                                    int NUM_sep_edges,
                                    int NUM_eles,
                                    int NUM_glob_interior_edges,
                                    T P_in,
                                    T P_out);
}; // namespace cuDFNsys