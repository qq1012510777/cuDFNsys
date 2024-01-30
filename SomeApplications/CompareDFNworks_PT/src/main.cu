#include "cuDFNsys.cuh"
#include <fstream>
#include <iostream>
#include <limits.h>
#include <numeric>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>

int main(int argc, char *argv[])
{
    time_t t;
    time(&t);
    string csvname("example_FracPara_PercolativeFractures");
    double minMeshSize = 0.1;
    double maxMeshSize = 0.17;

    cuDFNsys::DFN<double> my_dfn;
    cuDFNsys::MeshDFN<double> meshGen;
    cuDFNsys::FlowDFN<double> flowDFN;
    cuDFNsys::PTDFN<double> particleTracking;

    // my_dfn.RandomSeed = (unsigned long)t;
    my_dfn.LoadDFNFromCSV(csvname);
    my_dfn.IdentifyIntersectionsClusters(true);
    my_dfn.Visualization("DFN_VISUAL", "DFN_VISUAL", "DFN_VISUAL", true, false,
                         true, true);

    meshGen.MinElementSize = minMeshSize;
    meshGen.MaxElementSize = maxMeshSize;

    meshGen.MeshGeneration(my_dfn);
    meshGen.Visualization(my_dfn, "DFN_MESH_VISUAL", "DFN_MESH_VISUAL",
                          "DFN_MESH_VISUAL", true, true);

    //cout << "\n\n--------------------\nDFN gen and mesh Elapse time: "
    //     << cuDFNsys::CPUSecond() - iStart_DFN << " seconds\n\n";

    flowDFN.MuOverRhoG = 1;
    flowDFN.InletHead = 30;
    flowDFN.OutletHead = 0;

    flowDFN.FlowSimulation(my_dfn, meshGen);
    flowDFN.Visualization(my_dfn, meshGen, "DFN_FLOW_VISUAL", "DFN_FLOW_VISUAL",
                          "DFN_FLOW_VISUAL");

    //cout << "\n\n--------------------\nDFN flow Elapse time: "
    //     << cuDFNsys::CPUSecond() - iStart_DFN << " seconds\n\n";

    particleTracking.NumParticles = 24000;
    particleTracking.NumTimeSteps = 1000;
    particleTracking.MolecularDiffusion = 0;
    particleTracking.DeltaT =
        0.2 * pow(meshGen.MeanGridSize, 0.5) / flowDFN.MaxVelocity;
    particleTracking.InjectionMethod = "Resident";
    particleTracking.OutputAllPTInformationOrFPTCurve = false;
    particleTracking.SpacingOfControlPlanes = 100000;
    particleTracking.IfOutputVarianceOfDisplacementsEachStep = true;
    particleTracking.IfInjectAtCustomedPlane = true;
    particleTracking.CustomedPlaneInjection = 14;
    particleTracking.IfUseFluxWeightedOrEqualProbableMixingIntersection = true;

    double iStart_DFN = cuDFNsys::CPUSecond();
    particleTracking.ParticleTracking(my_dfn, meshGen, flowDFN);
    cout << "\n\n--------------------\nDFN PT Elapse time: "
         << cuDFNsys::CPUSecond() - iStart_DFN << " seconds\n\n";

    particleTracking.Visualization(
        my_dfn, meshGen, flowDFN, "DFN_DISPERSION_VISUAL",
        "DFN_DISPERSION_VISUAL", "DFN_FLOW_VISUAL_I");

    return 0;
}