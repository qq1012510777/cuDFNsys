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

// ====================================================
// NAME:        A quickstart example to generate DFNs
// DESCRIPTION: Call cuDFNsys functions to do simulation.
// AUTHOR:      Tingchang YIN
// DATE:        13/10/2023
// ====================================================

#include "cuDFNsys.cuh"
#include <fstream>
#include <iostream>
#include <limits.h>
#include <unistd.h>

int main(int argc, char *argv[])
{
    time_t t;
    time(&t);

    cuDFNsys::DFN<double> my_dfn;

    my_dfn.NumFractures = {70, 80};
    my_dfn.Kappa = {20, 10};
    my_dfn.MeanOrientationOfFisherDistribution = {make_double3(0., 0., 1.),
                                                  make_double3(1., 0., 0.)};
    my_dfn.DomainSizeX = 30;
    my_dfn.DomainDimensionRatio = make_double3(1., 1., 2.);
    my_dfn.Beta = {0.2, 0.3};
    my_dfn.Gamma = {2.0e-5, 3.0e-6};
    my_dfn.ModeOfSizeDistribution = {0, 1};
    my_dfn.SizeDistributionParameters = {make_double4(1.5, 1., 15., 0.),
                                         make_double4(8.5, 5.5, 1., 15.)};
    my_dfn.PercoDir = 2;
    my_dfn.RandomSeed = (unsigned long)t;
    my_dfn.FractureGeneration();
    my_dfn.IdentifyIntersectionsClusters(true);
    my_dfn.Visualization("DFN_VISUAL", "DFN_VISUAL", "DFN_VISUAL", true, true,
                         true, true);

    cuDFNsys::MeshDFN<double> meshGen;
    meshGen.MinElementSize = 1;
    meshGen.MaxElementSize = 3;
    meshGen.MeshGeneration(my_dfn);
    meshGen.Visualization(my_dfn, "DFN_MESH_VISUAL", "DFN_MESH_VISUAL",
                          "DFN_MESH_VISUAL", true, true);

    cuDFNsys::FlowDFN<double> flowDFN;
    flowDFN.InletHead = 60;
    flowDFN.OutletHead = 0;
    flowDFN.FlowSimulation(my_dfn, meshGen);
    flowDFN.Visualization(my_dfn, meshGen, "DFN_FLOW_VISUAL", "DFN_FLOW_VISUAL",
                          "DFN_FLOW_VISUAL");

    cuDFNsys::PTDFN<double> particleTracking;
    particleTracking.NumParticles = 20000;
    particleTracking.NumTimeSteps = 200;
    particleTracking.PecletNumber = 300;
    particleTracking.LengthScalePe = 30;
    particleTracking.VelocityScalePe = flowDFN.MeanVelocity;
    particleTracking.MolecularDiffusion = particleTracking.LengthScalePe /
                                          particleTracking.PecletNumber *
                                          particleTracking.VelocityScalePe;
    particleTracking.FactorTimeScaleCrossElement = 2;
    particleTracking.TimeScaleCrossElement =
        pow(meshGen.MeanGridSize, 0.5) / flowDFN.MaxVelocity;
    particleTracking.DeltaT = particleTracking.TimeScaleCrossElement /
                              particleTracking.FactorTimeScaleCrossElement;
    particleTracking.FluexWeightedOrUniformInjection = true;
    particleTracking.OutputAllPTInformationOrFPTCurve = true;
    particleTracking.SpacingOfControlPlanes = 30;
    particleTracking.IfOutputVarianceOfDisplacementsEachStep = true;
    particleTracking.IfInjectAtCustomedPlane = true;
    particleTracking.CustomedPlaneInjection = 23;
    particleTracking.IfUseFluxWeightedOrEqualProbableMixingIntersection = true;

    particleTracking.ParticleTracking(my_dfn, meshGen, flowDFN);
    particleTracking.Visualization(my_dfn, meshGen, flowDFN,
                                   "DFN_DISPERSION_VISUAL",
                                   "DFN_DISPERSION_VISUAL", "DFN_FLOW_VISUAL");
    return 0;
};
