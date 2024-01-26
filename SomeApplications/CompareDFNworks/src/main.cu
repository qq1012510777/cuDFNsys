#include "cuDFNsys.cuh"
#include <string>

int main(int argc, char *argv[])
{
    time_t t;
    time(&t);
    cuDFNsys::DFN<double> my_dfn;
    cuDFNsys::MeshDFN<double> meshGen;
    cuDFNsys::FlowDFN<double> flowDFN;

    string csvname(argv[1]);

    double iStart_DFN = cuDFNsys::CPUSecond();
    my_dfn.LoadDFNFromCSV(csvname);
    my_dfn.IdentifyIntersectionsClusters(true);
    my_dfn.Visualization("DFN_VISUAL", "DFN_VISUAL", "DFN_VISUAL", true, false,
                         true, true);

    meshGen.MinElementSize = atof(argv[2]);
    meshGen.MaxElementSize = atof(argv[3]);
    meshGen.MeshGeneration(my_dfn);
    meshGen.Visualization(my_dfn, "DFN_MESH_VISUAL", "DFN_MESH_VISUAL",
                          "DFN_MESH_VISUAL", true, true);

    cout << "\n\n--------------------\nDFN gen and mesh Elapse time: "
         << cuDFNsys::CPUSecond() - iStart_DFN << " seconds\n\n";

    flowDFN.MuOverRhoG = 1;
    flowDFN.InletHead = 30;
    flowDFN.OutletHead = 0;

    iStart_DFN = cuDFNsys::CPUSecond();
    flowDFN.FlowSimulation(my_dfn, meshGen);
    flowDFN.Visualization(my_dfn, meshGen, "DFN_FLOW_VISUAL", "DFN_FLOW_VISUAL",
                          "DFN_FLOW_VISUAL");

    cout << "\n\n--------------------\nDFN flow Elapse time: "
         << cuDFNsys::CPUSecond() - iStart_DFN << " seconds\n\n";
    return 0;
}