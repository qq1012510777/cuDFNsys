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
#include "CheckIfReachControlPlanesKernel.cuh"
#include "EdgeToEle.cuh"
#include "IsNotEqual.cuh"
#include "NeighborEle.cuh"
#include "OutputObjectData/OutputObjectData.cuh"
#include "Particle.cuh"
#include "ParticleMovementOneTimeStepCPU.cuh"
#include "ParticleMovementOneTimeStepGPUKernel.cuh"
#include "PredicateNumOfReachedOutletParticles.cuh"
#include "ThrustStatistics/ThrustStatistics.cuh"
#include "ToStringWithWidth/ToStringWithWidth.cuh"
#include "Transform2DTo3DKernel.cuh"
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <thrust/count.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/remove.h>
#include <thrust/sort.h>
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
        string RecordMode = "OutputAll";

        thrust::host_vector<T> ControlPlanes;

        bool IfOutputMSD = true;

        thrust::host_vector<int2> CorrespondingEleLocalEdge;

        bool IfPeriodic = false;

        uint TimeIntervalOutPTInformation;

        //------------record the travel time reaching control planes
        thrust::host_vector<uint> TimeReachControlPlanes;

    private:
        string ParticlePosition = "ParticlePositionResult/ParticlePosition";
        string DispersionInfo = "ParticlePositionResult/DispersionInfo";

    public:
        ParticleTransport(){};
        ///   ParticleTransport(const int &NumOfParticles,
        ///                     const int &NumTimeStep,
        ///                     T delta_T_,
        ///                     T Dispersion_local,
        ///                     thrust::host_vector<cuDFNsys::Fracture<T>> Fracs,
        ///                     cuDFNsys::Mesh<T> mesh,
        ///                     const cuDFNsys::MHFEM<T> &fem,
        ///                     uint Dir_flow,
        ///                     T outletcoordinate,
        ///                     const string &Particle_mode,
        ///                     const string &Injection_mode,
        ///                     bool if_cpu = false,
        ///                     int Nproc = 10,
        ///                     bool record_time = false, // record the run time of each step
        ///                     string recordMode = "OutputAll");

        ParticleTransport(const int &NumTimeStep,
                          thrust::host_vector<cuDFNsys::Fracture<T>> Fracs,
                          cuDFNsys::Mesh<T> mesh, const cuDFNsys::MHFEM<T> &fem,
                          uint Dir_flow, T outletcoordinate,
                          int NumOfParticles_ii = 0, T delta_T_ii = 0,
                          T Diffusion_local_ii = 0, // molecular diffusion
                          string Particle_mode_ii = "Particle_tracking",
                          string Injection_mode_ii = "Flux-weighted",
                          string recordMode = "OutputAll", bool if_cpu = false,
                          int Nproc = 10,
                          bool record_time = false, // record running time
                          T SpacingOfControlPlanes = 10,
                          bool IfOutputMSD = true,
                          bool IfInitCenterDomain = false, T InjectionPlane = 0,
                          bool If_completeMixing_fluxWeighted = true,
                          bool IfPeriodic_ = false,
                          uint TimeIntervalOutPTInformation_s = 100);

        void ParticleMovement(const int &init_NO_STEP, const int &NumTimeStep,
                              T delta_T_, T Dispersion_local,
                              const string &Particle_mode,
                              const string &Injection_mode,
                              thrust::host_vector<cuDFNsys::Fracture<T>> Fracs,
                              cuDFNsys::Mesh<T> mesh,
                              const cuDFNsys::MHFEM<T> &fem, T outletcoordinate,
                              bool If_completeMixing_fluxWeighted = true);

        void ParticleMovementCPU(
            const int &init_NO_STEP, const int &NumTimeStep, T delta_T_,
            T Dispersion_local, const string &Particle_mode,
            const string &Injection_mode,
            thrust::host_vector<cuDFNsys::Fracture<T>> Fracs,
            cuDFNsys::Mesh<T> mesh, const cuDFNsys::MHFEM<T> &fem,
            T outletcoordinate, int Nproc = 10);

        void OutputParticleInfoStepByStep(
            const uint &StepNO, const T delta_T, const T Dispersion_local,
            const string &Particle_mode, const string &Injection_mode,
            thrust::host_vector<cuDFNsys::Fracture<T>> Fracs,
            cuDFNsys::Mesh<T> mesh);

        void OutputMSD(const uint &StepNO,
                       const thrust::host_vector<cuDFNsys::Fracture<T>> &Fracs,
                       const cuDFNsys::Mesh<T> &mesh,
                       const T &HalfDomainSize_PercoDirection);

        thrust::host_vector<T> Get3DParticlePositions(
            const thrust::host_vector<cuDFNsys::Fracture<T>> &Fracs,
            const cuDFNsys::Mesh<T> &mesh, T L_percoDir);

        void IfReachControlPlane(
            const uint &StepNo, const uint &PercoDir,
            const std::vector<T> &ControlPlane,
            const thrust::host_vector<cuDFNsys::Fracture<T>> &Fracs,
            const cuDFNsys::Mesh<T> &mesh, const T &HalfDomainSize_PercoDir);

        void MatlabPlot(const string &mat_key, const string &command_key,
                        const cuDFNsys::Mesh<T> &mesh,
                        const cuDFNsys::MHFEM<T> &fem, const T &L,
                        double3 DomainDimensionRatio = make_double3(1, 1, 1),
                        bool if_python_visualization = false,
                        string PythonName_Without_suffix = "ParticleMovement");

    private:
        void IdentifyEdgesSharedEle(cuDFNsys::Mesh<T> mesh);

        void InitilizeParticles(const int &NumOfParticles,
                                cuDFNsys::Mesh<T> mesh,
                                const cuDFNsys::MHFEM<T> &fem,
                                const string &Injection_mode,
                                bool IfInitCenterDomain = false,
                                T InjectionPlane = 0);

        void IdentifyNeighborElements(cuDFNsys::Mesh<T> mesh);

        void
        IdentifyInletOutletCorrespondingElementEdge(cuDFNsys::Mesh<T> mesh);
    };
}; // namespace cuDFNsys