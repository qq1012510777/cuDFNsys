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
#include "../DataTypeSelector/DataTypeSelector.cuh"
#include "../GlobalDef/GlobalDef.cuh"
#include "../HDF5API/HDF5API.cuh"
#include "../Mesh/Mesh.cuh"
#include "AssembleOnGPUKernel.cuh"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/UmfPackSupport"
#include "ReconstructVelocityGrid.cuh"
#include "Triplet.cuh"

using namespace Eigen;

namespace cuDFNsys
{
    template <typename T>
    class MHFEM
    {
    public:
        T QIn = 0;
        T QOut = 0;
        T InletLength = 0;
        T OutletLength = 0;
        T QError = 0;
        T Permeability = 0;
        Eigen::MatrixXd PressureInteriorEdge;
        Eigen::MatrixXd PressureEles;
        Eigen::MatrixXd VelocityNormalScalarSepEdges;
        T InletP = 100;
        T OutletP = 20;
        double TripletTime = 0;
        double MaxVelocity;
        double MeanVelocity;

    public:
        T MuOverRhoG = 1;
        bool IfPeriodic = false;
        T ConsTq = 1e-15;

    public:
        int Dir = 2;

    public:
        MHFEM(){};

        MHFEM(const cuDFNsys::Mesh<T> &mesh,
              const thrust::host_vector<cuDFNsys::Fracture<T>> &Fracs,
              const T &inlet_p_, const T &outlet_p_, const int &dir_,
              const T &L, double3 DomainDimensionRatio = make_double3(1, 1, 1),
              bool if_CPU = false, int Nproc = 10, bool if_periodic = false,
              T muOverRhoG = 1, T constant_qq = 1e-5);

        double2 MatlabPlot(const string &mat_key, const string &command_key,
                           thrust::host_vector<cuDFNsys::Fracture<T>> Fracs,
                           const cuDFNsys::Mesh<T> &mesh, const T &L,
                           bool if_python_visualization = false,
                           string PythonName_Without_suffix = "DFN_mhfem_py",
                           double3 DomainDimensionRatio = make_double3(1, 1,
                                                                       1));

    private:
        void
        Implementation(const cuDFNsys::Mesh<T> &mesh,
                       const thrust::host_vector<cuDFNsys::Fracture<T>> &Fracs,
                       bool if_CPU = false, int Nproc = 10);

    private:
        pair<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>>
        AssembleOnGPU(const cuDFNsys::Mesh<T> &mesh,
                      const thrust::host_vector<cuDFNsys::Fracture<T>> &Fracs,
                      T P_in, T P_out);

        pair<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>>
        AssembleOnCPU(const cuDFNsys::Mesh<T> &mesh,
                      const thrust::host_vector<cuDFNsys::Fracture<T>> &Fracs,
                      T P_in, T P_out, int Nproc = 10);
    };

}; // namespace cuDFNsys