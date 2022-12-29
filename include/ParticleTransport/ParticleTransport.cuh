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
// NAME:              ParticleTransport.cuh
//
// PURPOSE:           perform particle tracking with random walk method
//
// FUNCTIONS/OBJECTS: ParticleTransport
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "EdgeToEle.cuh"
#include "NeighborEle.cuh"
#include "OutputObjectData/OutputObjectData.cuh"
#include "Particle.cuh"
#include "ParticleMovementOneTimeStepCPU.cuh"
#include "ParticleMovementOneTimeStepGPUKernel.cuh"
#include "PredicateNumOfReachedOutletParticles.cuh"
#include "ToStringWithWidth/ToStringWithWidth.cuh"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

namespace cuDFNsys
{
template <typename T>
class ParticleTransport
{
public:
    int NumParticles = 0;
    thrust::host_vector<cuDFNsys::Particle<T>> ParticlePlumes;
    thrust::host_vector<cuDFNsys::EdgeToEle> EdgesSharedEle;
    thrust::host_vector<cuDFNsys::NeighborEle> NeighborEleOfOneEle;
    uint Dir = 2;
    uint SizeOfDataBlock = 2000; // how many steps store in a h5 file.
    uint BlockNOPresent = 0;
    vector<double> RunTimeEveryStep;
    bool IfRecordTime = false;

private:
    string ParticlePosition = "ParticlePositionResult/ParticlePosition";
    string DispersionInfo = "ParticlePositionResult/DispersionInfo";

public:
    ParticleTransport(const int &NumOfParticles,
                      const int &NumTimeStep,
                      T delta_T_,
                      T Dispersion_local,
                      thrust::host_vector<cuDFNsys::Fracture<T>> Fracs,
                      cuDFNsys::Mesh<T> mesh,
                      const cuDFNsys::MHFEM<T> &fem,
                      uint Dir_flow,
                      T outletcoordinate,
                      const string &Particle_mode,
                      const string &Injection_mode,
                      bool if_cpu = false, int Nproc = 10, bool record_time = false);

    ParticleTransport(const int &NumTimeStep,
                      thrust::host_vector<cuDFNsys::Fracture<T>> Fracs,
                      cuDFNsys::Mesh<T> mesh,
                      const cuDFNsys::MHFEM<T> &fem,
                      uint Dir_flow,
                      T outletcoordinate,
                      int NumOfParticles_ii = 0,
                      T delta_T_ii = 0,
                      T Dispersion_local_ii = 0,
                      string Particle_mode_ii = "Particle_tracking",
                      string Injection_mode_ii = "Flux-weighted",
                      bool if_cpu = false, int Nproc = 10, bool record_time = false);

    void ParticleMovement(const int &init_NO_STEP,
                          const int &NumTimeStep,
                          T delta_T_,
                          T Dispersion_local,
                          const string &Particle_mode,
                          const string &Injection_mode,
                          thrust::host_vector<cuDFNsys::Fracture<T>> Fracs,
                          cuDFNsys::Mesh<T> mesh,
                          const cuDFNsys::MHFEM<T> &fem,
                          T outletcoordinate);

    void ParticleMovementCPU(const int &init_NO_STEP,
                             const int &NumTimeStep,
                             T delta_T_,
                             T Dispersion_local,
                             const string &Particle_mode,
                             const string &Injection_mode,
                             thrust::host_vector<cuDFNsys::Fracture<T>> Fracs,
                             cuDFNsys::Mesh<T> mesh,
                             const cuDFNsys::MHFEM<T> &fem,
                             T outletcoordinate,
                             int Nproc = 10);

    void OutputParticleInfoStepByStep(const uint &StepNO,
                                      const T delta_T,
                                      const T Dispersion_local,
                                      const string &Particle_mode,
                                      const string &Injection_mode,
                                      thrust::host_vector<cuDFNsys::Fracture<T>> Fracs,
                                      cuDFNsys::Mesh<T> mesh);

    void MatlabPlot(const string &mat_key,
                    const string &command_key,
                    const cuDFNsys::Mesh<T> &mesh,
                    const cuDFNsys::MHFEM<T> &fem,
                    const T &L,
                    bool if_python_visualization = false,
                    string PythonName_Without_suffix = "ParticleMovement");

private:
    void IdentifyEdgesSharedEle(cuDFNsys::Mesh<T> mesh);

    void InitilizeParticles(const int &NumOfParticles,
                            cuDFNsys::Mesh<T> mesh,
                            const cuDFNsys::MHFEM<T> &fem,
                            const string &Injection_mode);

    void IdentifyNeighborElements(cuDFNsys::Mesh<T> mesh);

};
}; // namespace cuDFNsys