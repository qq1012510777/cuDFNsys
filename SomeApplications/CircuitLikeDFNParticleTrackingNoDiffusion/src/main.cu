// more details of this code is shown in the manual
// more details of this code is shown in the manual
// more details of this code is shown in the manual
#include "cuDFNsys.cuh"
int main()
{
    time_t t;
    time(&t);

    cuDFNsys::DFN<double> my_dfn;
    cuDFNsys::MeshDFN<double> meshGen;
    cuDFNsys::FlowDFN<double> flowDFN;
    cuDFNsys::PTDFN<double> particleTracking;

    my_dfn.RandomSeed = (unsigned long)t;
    my_dfn.LoadDFNFromCSV("InputParametersForDeterministicDFN");

    my_dfn.IdentifyIntersectionsClusters(true);
    my_dfn.Visualization("DFN_VISUAL", "DFN_VISUAL", "DFN_VISUAL", true, true,
                         true, true);

    meshGen.LoadParametersFromCSV("Mesh_parameters");
    meshGen.MeshGeneration(my_dfn);
    meshGen.Visualization(my_dfn, "DFN_MESH_VISUAL", "DFN_MESH_VISUAL",
                          "DFN_MESH_VISUAL", true, true);

    flowDFN.LoadParametersFromCSV("Flow_parameters");
    flowDFN.FlowSimulation(my_dfn, meshGen);
    flowDFN.Visualization(my_dfn, meshGen, "DFN_FLOW_VISUAL", "DFN_FLOW_VISUAL",
                          "DFN_FLOW_VISUAL");

    int u = 0;
    cout << "If run PT?\n";
    cin >> u;

    if (u != 0)
    {
        particleTracking.LoadParametersFromCSV("PT_parameters");
        particleTracking.ParticleTracking(my_dfn, meshGen, flowDFN);
        particleTracking.Visualization(
            my_dfn, meshGen, flowDFN, "DFN_DISPERSION_VISUAL",
            "DFN_DISPERSION_VISUAL", "DFN_FLOW_VISUAL");
    }
    my_dfn.StoreInH5("Class_DFN");
    meshGen.StoreInH5("Class_MESH");
    flowDFN.StoreInH5("Class_FLOW");

    // cout << "*** Right now the particle data are two-dimensional. Use "
    //      << "the executable `Transform2DH5ParticleDataTo3D` to "
    //      << "transform "
    //      << "them to 3D!  ***" << endl;
    // cout << "*** Just run: ./Transform2DH5ParticleDataTo3D 0 "
    //      << "DFN_MESH_VISUAL.h5 ***" << std::endl;
    return 0;
};
