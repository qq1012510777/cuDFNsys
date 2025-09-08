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
// NAME:              cuMechsysyDFN.cuh
//
// PURPOSE:           The API of cuMechsysyDFN
//
// FUNCTIONS/OBJECTS: N/A
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////

#pragma once
#include "./CPUSecond/CPUSecond.cuh"
#include "./DataTypeSelector/DataTypeSelector.cuh"
#include "./Exceptions/Exceptions.cuh"
#include "./Fractures/Fracture.cuh"
#include "./Fractures/Fractures.cuh"
#include "./Fractures/GetAllPercolatingFractures.cuh"
#include "./Fractures/IdentifyIntersection.cuh"
#include "./Fractures/IdentifyIntersectionKernel.cuh"
#include "./Fractures/IdentifyPercolationCluster.cuh"
#include "./Fractures/Intersection.cuh"
#include "./Fractures/MatlabPlotDFN.cuh"
#include "./Fractures/RemoveDeadEndFrac.cuh"
#include "./GPUErrCheck/GPUErrCheck.cuh"
#include "./Geometry/Geometry.cuh"
#include "./GetStatistics/GetStatistics.cuh"
#include "./GlobalDef/GlobalDef.cuh"
#include "./Graph/Graph.cuh"
#include "./HDF5API/HDF5API.cuh"
#include "./MHFEM/AssembleOnGPUKernel.cuh"
#include "./MHFEM/MHFEM.cuh"
#include "./MHFEM/ReconstructVelocityGrid.cuh"
#include "./MHFEM/StimaA.cuh"
#include "./MHFEM/Triplet.cuh"
// #include "./MatlabAPI/MatlabAPI.cuh"
#include "./MatrixManipulation/MatrixManipulation.cuh"
#include "./Mesh/Mesh.cuh"
#include "./Quaternion/Quaternion.cuh"
#include "./RandomFunction/RandomFunction.cuh"
#include "./ToStringWithWidth/ToStringWithWidth.cuh"
#include "./Warmup/Warmup.cuh"

#include "./ParticleTransport/EdgeToEle.cuh"
#include "./ParticleTransport/IdentifyParticleCrossesWhichEdge.cuh"
#include "./ParticleTransport/IfParticleOnBoundOfElement.cuh"
#include "./ParticleTransport/IfParticlePositionNearOneVertexOfElement.cuh"
#include "./ParticleTransport/NeighborEle.cuh"
#include "./ParticleTransport/Particle.cuh"
#include "./ParticleTransport/ParticleMovementOneTimeStepGPUKernel.cuh"
#include "./ParticleTransport/ParticleReflection.cuh"
#include "./ParticleTransport/ParticleTransport.cuh"
#include "./ParticleTransport/Roate2DPositionTo3D.cuh"
#include "./ParticleTransport/WhichElementToGo.cuh"

#include "./InputObjectData/InputObjectData.cuh"
#include "./OutputObjectData/OutputObjectData.cuh"

#include <thrust/count.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

// #include "./CPUSecond/CPUSecond.cuh"
// #include "./Exceptions/Exceptions.cuh"
// #include "./Fractures/Fractures.cuh"
// #include "./Fractures/GetAllPercolatingFractures.cuh"
// #include "./Fractures/IdentifyIntersection.cuh"
// #include "./Fractures/IdentifyPercolationCluster.cuh"
// #include "./Fractures/MatlabPlotDFN.cuh"
// #include "./Fractures/RemoveDeadEndFrac.cuh"
// #include "./GPUErrCheck/GPUErrCheck.cuh"
// #include "./GetStatistics/GetStatistics.cuh"
// #include "./GlobalDef/GlobalDef.cuh"
// #include "./Graph/Graph.cuh"
// #include "./HDF5API/HDF5API.cuh"
// #include "./MHFEM/MHFEM.cuh"
// #include "./MatlabAPI/MatlabAPI.cuh"
// #include "./Mesh/Mesh.cuh"
// #include "./ParticleTransport/EdgeToEle.cuh"
// #include "./ParticleTransport/Particle.cuh"
// #include "./ParticleTransport/ParticleTransport.cuh"
// #include "./Quaternion/Quaternion.cuh"
// #include "./ToStringWithWidth/ToStringWithWidth.cuh"
// #include "./Warmup/Warmup.cuh"

namespace cuDFNsys
{
    template <typename T>
    class DFN
    {
    public:
        std::vector<int> NumFractures;
        std::vector<T> Kappa;
        std::vector<cuDFNsys::Vector3<T>> MeanOrientationOfFisherDistribution;
        T DomainSizeX;
        double3 DomainDimensionRatio;
        std::vector<T> Beta;
        std::vector<T> Gamma;
        std::vector<int> ModeOfSizeDistribution;
        std::vector<cuDFNsys::Vector4<T>> SizeDistributionParameters;
        int PercoDir;
        int NumFracturesTotal;
        unsigned long RandomSeed;

        thrust::host_vector<cuDFNsys::Fracture<T>> FracturesHost;
        thrust::device_vector<cuDFNsys::Fracture<T>> FracturesDevice;
        cuDFNsys::Fracture<T> *FracturesDevicePtr;
        std::map<pair<size_t, size_t>,
                 pair<cuDFNsys::Vector3<T>, cuDFNsys::Vector3<T>>>
            IntersectionMap;
        std::vector<std::vector<size_t>> ListClusters;
        std::vector<size_t> PercolationCluster;

        bool IfPseudo3D = false;

    private:
        bool IfPeriodic = false;

    public:
        DFN() {};
        void FractureGeneration();
        void IdentifyIntersectionsClusters(const bool &IfTruncatedFractures);
        void Visualization(const string &MatlabScriptName,
                           const string &PythonScriptName,
                           const string &HDF5FileName,
                           const bool &IfShowTruncatedFractures,
                           const bool &IfShowIntersections,
                           const bool &IfHightlighAllClusters,
                           const bool &IfShowOrientationDistribution);
        void StoreInH5(const string &ClassNameH5);
        void LoadClassFromH5(const string &ClassNameH5);
        void SpatialPeriodicity();
        void LoadDFNFromCSV(const string &xlsxNameWithoutSuffix);
        bool CheckIfPeriodic() { return this->IfPeriodic; };
        void ChangeDomainSize(const T &L, const double3 &DomainDimensionRatioTT);
        void ChangePercolationDirectionIdentifyPercolationCluster(const int &PerDir);
    };
}; // namespace cuDFNsys

namespace cuDFNsys
{
    template <typename T>
    class MeshDFN
    {
    public:
        cuDFNsys::Mesh<T> MeshData;
        T MinElementSize;
        T MaxElementSize;
        T MeanGridSize;
        std::vector<size_t> FracsPercol;

    public:
        MeshDFN() {};
        void MeshGeneration(cuDFNsys::DFN<T> &my_dfn);
        void Visualization(cuDFNsys::DFN<T> my_dfn,
                           const string &MatlabScriptName,
                           const string &PythonScriptName,
                           const string &HDF5FileName,
                           const bool &IfCheck2DCoordinatesOfMesh,
                           const bool &IfCheckEdgeAttributes);
        void StoreInH5(const string &ClassNameH5);
        void LoadClassFromH5(const string &ClassNameH5);
        void LoadParametersFromCSV(const string &CSVName);
        void ChangePecolationDirectionAndRenumberingEdge(const int PercoDir, const T L, double3 DomainDimensionRatio = make_double3(1, 1, 1));
    };
}; // namespace cuDFNsys

namespace cuDFNsys
{
    template <typename T>
    class FlowDFN
    {
    public:
        T InletHead;
        T OutletHead;
        T MaxVelocity;
        T MeanVelocity;
        cuDFNsys::MHFEM<T> FlowData;
        T MuOverRhoG = 1;
        T ConsTq = 1e-15;

    public:
        bool IfPeriodic = false;

    public:
        FlowDFN() {};
        void FlowSimulation(cuDFNsys::DFN<T> my_dfn,
                            cuDFNsys::MeshDFN<T> my_mesh);
        void Visualization(cuDFNsys::DFN<T> my_dfn,
                           cuDFNsys::MeshDFN<T> my_mesh,
                           const string &MatlabScriptName,
                           const string &PythonScriptName,
                           const string &HDF5FileName);
        void StoreInH5(const string &ClassNameH5);
        void LoadClassFromH5(const string &ClassNameH5);
        void LoadParametersFromCSV(const string &CSVName);
    };
}; // namespace cuDFNsys

namespace cuDFNsys
{
    template <typename T>
    class PTDFN
    {
    public:
        int NumParticles;
        int NumTimeSteps;
        T PecletNumber;
        T LengthScalePe;
        T VelocityScalePe;
        T MolecularDiffusion;
        T FactorTimeScaleCrossElement;
        T DeltaT;
        T TimeScaleCrossElement;
        string InjectionMethod;
        bool IfUseFluxWeightedOrEqualProbableMixingIntersection;
        T SpacingOfControlPlanes;
        bool IfOutputVarianceOfDisplacementsEachStep;
        bool IfInjectAtCustomedPlane;
        T CustomedPlaneInjection;
        bool OutputAllPTInformationOrFPTCurve;
        cuDFNsys::ParticleTransport<T> PTData;
        bool IfPeriodic = false;
        uint TimeIntervalOutPTInformation = 100;
        bool IfOutputAllParticleAccumulativeDisplacement = false;
        size_t IfPureDiffusion = 0; 
        size_t IfDiscontinueAfterFirstAbsorption = 0;
        size_t IfReflectionAtInlet = 0;
    public:
        PTDFN() {};
        void ParticleTracking(cuDFNsys::DFN<T> my_dfn,
                              cuDFNsys::MeshDFN<T> my_mesh,
                              cuDFNsys::FlowDFN<T> my_flow);
        void Visualization(cuDFNsys::DFN<T> my_dfn,
                           cuDFNsys::MeshDFN<T> my_mesh,
                           cuDFNsys::FlowDFN<T> my_flow,
                           const string &MatlabScriptName,
                           const string &PythonScriptName,
                           const string &HDF5FileNameOfFlowDFN);
        void LoadParametersFromCSV(const string &CSVName);
    };
}; // namespace cuDFNsys