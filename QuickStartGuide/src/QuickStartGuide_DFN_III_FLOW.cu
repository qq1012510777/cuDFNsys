#include "cuDFNsys.cuh"
int main()
{
    cuDFNsys::DFN<double> my_dfn;
    my_dfn.LoadClassFromH5("Class_DFN");
    cuDFNsys::MeshDFN<double> meshGen;
    meshGen.LoadClassFromH5("Class_MESH");
    meshGen.Visualization(my_dfn, "DFN_MESH_VISUAL_II", "DFN_MESH_VISUAL_II",
                          "DFN_MESH_VISUAL_II", true, true);
    cuDFNsys::FlowDFN<double> flowDFN;
    flowDFN.InletHead = 60;
    flowDFN.OutletHead = 0;
    flowDFN.FlowSimulation(my_dfn, meshGen);
    flowDFN.Visualization(my_dfn, meshGen, "DFN_FLOW_VISUAL_I",
                          "DFN_FLOW_VISUAL_I", "DFN_FLOW_VISUAL_I");
    flowDFN.StoreInH5("Class_FLOW");
    return 0;
};