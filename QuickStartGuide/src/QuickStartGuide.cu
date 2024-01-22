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

    cout << "Please input the minimum grid size you wanted:\n";
    cin >> meshGen.MinElementSize;
    cout << "Please input the maximum grid size you wanted:\n";
    cin >> meshGen.MaxElementSize;

    cout << "Please input the inlet head value (in meters) you wanted:\n";
    cin >> flowDFN.InletHead;
    cout << "Please input the outlet head value (in meters) you wanted:\n";
    cin >> flowDFN.OutletHead;

    my_dfn.RandomSeed = (unsigned long)t;
    my_dfn.LoadDFNFromCSV("InputParametersForStochasticDFN");

    my_dfn.IdentifyIntersectionsClusters(true);
    my_dfn.Visualization("DFN_VISUAL", "DFN_VISUAL", "DFN_VISUAL", true, true,
                         true, true);

    meshGen.MeshGeneration(my_dfn);
    meshGen.Visualization(my_dfn, "DFN_MESH_VISUAL", "DFN_MESH_VISUAL",
                          "DFN_MESH_VISUAL", true, true);

    flowDFN.FlowSimulation(my_dfn, meshGen);
    flowDFN.Visualization(my_dfn, meshGen, "DFN_FLOW_VISUAL", "DFN_FLOW_VISUAL",
                          "DFN_FLOW_VISUAL");

    particleTracking.LoadParametersFromCSV("PT_parameters");

    particleTracking.ParticleTracking(my_dfn, meshGen, flowDFN);
    particleTracking.Visualization(my_dfn, meshGen, flowDFN,
                                   "DFN_DISPERSION_VISUAL",
                                   "DFN_DISPERSION_VISUAL", "DFN_FLOW_VISUAL");

    my_dfn.StoreInH5("Class_DFN");
    meshGen.StoreInH5("Class_MESH");
    flowDFN.StoreInH5("Class_FLOW");

    cout << "*** Right now the particle data are two-dimensional. Use "
         << "the executable `Transform2DH5ParticleDataTo3D` to "
         << "transform "
         << "them to 3D!  ***" << endl;
    cout << "*** Just run: ./Transform2DH5ParticleDataTo3D 0 "
         << "DFN_MESH_VISUAL.h5 ***" << std::endl;
    return 0;
};
