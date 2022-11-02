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
//#include "./MatlabAPI/MatlabAPI.cuh"
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

#include "./OutputObjectData/OutputObjectData.cuh"

#include "./InputObjectData/InputObjectData.cuh"

//#include "./CPUSecond/CPUSecond.cuh"
//#include "./Exceptions/Exceptions.cuh"
//#include "./Fractures/Fractures.cuh"
//#include "./Fractures/GetAllPercolatingFractures.cuh"
//#include "./Fractures/IdentifyIntersection.cuh"
//#include "./Fractures/IdentifyPercolationCluster.cuh"
//#include "./Fractures/MatlabPlotDFN.cuh"
//#include "./Fractures/RemoveDeadEndFrac.cuh"
//#include "./GPUErrCheck/GPUErrCheck.cuh"
//#include "./GetStatistics/GetStatistics.cuh"
//#include "./GlobalDef/GlobalDef.cuh"
//#include "./Graph/Graph.cuh"
//#include "./HDF5API/HDF5API.cuh"
//#include "./MHFEM/MHFEM.cuh"
//#include "./MatlabAPI/MatlabAPI.cuh"
//#include "./Mesh/Mesh.cuh"
//#include "./ParticleTransport/EdgeToEle.cuh"
//#include "./ParticleTransport/Particle.cuh"
//#include "./ParticleTransport/ParticleTransport.cuh"
//#include "./Quaternion/Quaternion.cuh"
//#include "./ToStringWithWidth/ToStringWithWidth.cuh"
//#include "./Warmup/Warmup.cuh"