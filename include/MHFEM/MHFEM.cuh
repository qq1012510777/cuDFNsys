///////////////////////////////////////////////////////////////////
// NAME:              MHFEM.cuh
//
// PURPOSE:           mixed hybrid finite element method
//
// FUNCTIONS/OBJECTS: MHFEM
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../GlobalDef/GlobalDef.cuh"
#include "../Mesh/Mesh.cuh"
#include "AssembleOnGPUKernel.cuh"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/UmfPackSupport"
#include "Triplet.cuh"

using namespace Eigen;

namespace cuDFNsys
{
class MHFEM
{
public:
    float QIn = 0;
    float QOut = 0;
    float InletLength = 0;
    float OutletLength = 0;
    float QError = 0;
    float Permeability = 0;
    Eigen::MatrixXd PressureInteriorEdge;
    Eigen::MatrixXd PressureEles;
    Eigen::MatrixXd VelocityNormalScalarSepEdges;

private:
    int Dir = 2;
    float InletP = 100;
    float OutletP = 20;

public:
    MHFEM(const cuDFNsys::Mesh &mesh,
          const thrust::host_vector<cuDFNsys::Fracture> &Fracs,
          const float &inlet_p_,
          const float &outlet_p_,
          const int &dir_,
          const float &L);

    void MatlabPlot(const string &mat_key,
                    const string &command_key,
                    const cuDFNsys::Mesh &mesh,
                    const float &L);

private:
    void Implementation(const cuDFNsys::Mesh &mesh,
                        const thrust::host_vector<cuDFNsys::Fracture> &Fracs);

private:
    pair<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>> AssembleOnGPU(const cuDFNsys::Mesh &mesh,
                                                                                 const thrust::host_vector<cuDFNsys::Fracture> &Fracs,
                                                                                 float P_in,
                                                                                 float P_out);
};

}; // namespace cuDFNsys