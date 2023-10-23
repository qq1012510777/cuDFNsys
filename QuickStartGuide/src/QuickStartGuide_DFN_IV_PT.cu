#include "cuDFNsys.cuh"
int main()
{
    cuDFNsys::DFN<double> my_dfn;
    my_dfn.LoadClassFromH5("Class_DFN");
    cuDFNsys::MeshDFN<double> meshGen;
    meshGen.LoadClassFromH5("Class_MESH");
    cuDFNsys::FlowDFN<double> flowDFN;
    flowDFN.LoadClassFromH5("Class_FLOW");
    flowDFN.Visualization(my_dfn, meshGen, "DFN_FLOW_VISUAL_II",
                          "DFN_FLOW_VISUAL_II", "DFN_FLOW_VISUAL_II");
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
    particleTracking.Visualization(
        my_dfn, meshGen, flowDFN, "DFN_DISPERSION_VISUAL",
        "DFN_DISPERSION_VISUAL", "DFN_FLOW_VISUAL_I");
    cout << "*** Right now the particle data are two-dimensional. Use "
         << "the executable `Transform2DH5ParticleDataTo3D` to "
         << "transform "
         << "them to 3D!  ***" << endl;
    cout << "*** Just run: ./Transform2DH5ParticleDataTo3D 0 "
         << "DFN_MESH_VISUAL_I.h5 ***" << std::endl;
    return 0;
};