#include "cuDFNsys.cuh"
int main()
{
    cuDFNsys::DFN<double> my_dfn;
    cuDFNsys::MeshDFN<double> meshGen;
    cuDFNsys::FlowDFN<double> flowDFN;
    cuDFNsys::PTDFN<double> particleTracking;

    my_dfn.LoadClassFromH5("Class_DFN");
    meshGen.LoadClassFromH5("Class_MESH");
    flowDFN.LoadClassFromH5("Class_FLOW");

    particleTracking.LoadParametersFromCSV("PT_parameters");
    particleTracking.ParticleTracking(my_dfn, meshGen, flowDFN);
    particleTracking.Visualization(my_dfn, meshGen, flowDFN,
                                   "DFN_DISPERSION_VISUAL",
                                   "DFN_DISPERSION_VISUAL", "DFN_FLOW_VISUAL");

    cout << "*** Right now the particle data are two-dimensional. Use "
         << "the executable `Transform2DH5ParticleDataTo3D` to "
         << "transform "
         << "them to 3D!  ***" << endl;
    cout << "*** Just run: ./Transform2DH5ParticleDataTo3D 0 "
         << "DFN_MESH_VISUAL.h5 ***" << std::endl;
    return 0;
};