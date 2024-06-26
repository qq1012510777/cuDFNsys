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

    double minMeshSize = 0.1;
    double maxMeshSize = 0.17;

    char cwd[PATH_MAX];
    if (getcwd(cwd, sizeof(cwd)) != NULL)
        printf("Current working dir: %s\n", cwd);
    else
        throw cuDFNsys::ExceptionsPause("getcwd() error");
    string curPath = cwd;

    int MCtime = 1;

    int Num = atoi(argv[1]);

    //string csvname("example_FracPara_PercolativeFractures");
    string csvname(argv[2]);

    std::vector<double> TimeDFN_GEN(MCtime), TimeDFN_MESH(MCtime),
        TimeDFN_FLOW(MCtime);

    for (int i = 0; i < MCtime; ++i)
    {
        chdir(curPath.c_str());
        string path2 = "DFN_" + cuDFNsys::ToStringWithWidth(Num, 3);
        string command1 = "mkdir -p " + path2;
        int resulty = system(command1.c_str());

        string command2 = curPath + "/" + path2;
        chdir(command2.c_str());

        cuDFNsys::DFN<double> my_dfn;
        cuDFNsys::MeshDFN<double> meshGen;
        cuDFNsys::FlowDFN<double> flowDFN;

        double iStart_DFN = cuDFNsys::CPUSecond();
        my_dfn.LoadDFNFromCSV((csvname));
        my_dfn.IdentifyIntersectionsClusters(true);
        my_dfn.Visualization("DFN_VISUAL", "DFN_VISUAL", "DFN_VISUAL", true,
                             false, true, true);
        TimeDFN_GEN[i] = cuDFNsys::CPUSecond() - iStart_DFN;

        meshGen.MinElementSize = minMeshSize;
        meshGen.MaxElementSize = maxMeshSize;

        iStart_DFN = cuDFNsys::CPUSecond();
        meshGen.MeshGeneration(my_dfn);
        meshGen.Visualization(my_dfn, "DFN_MESH_VISUAL", "DFN_MESH_VISUAL",
                              "DFN_MESH_VISUAL", true, true);

        //cout << "\n\n--------------------\nDFN gen and mesh Elapse time: "
        //     << cuDFNsys::CPUSecond() - iStart_DFN << " seconds\n\n";
        TimeDFN_MESH[i] = cuDFNsys::CPUSecond() - iStart_DFN;

        flowDFN.MuOverRhoG = 1;
        flowDFN.InletHead = 30;
        flowDFN.OutletHead = 0;

        iStart_DFN = cuDFNsys::CPUSecond();
        flowDFN.FlowSimulation(my_dfn, meshGen);
        flowDFN.Visualization(my_dfn, meshGen, "DFN_FLOW_VISUAL",
                              "DFN_FLOW_VISUAL", "DFN_FLOW_VISUAL");

        //cout << "\n\n--------------------\nDFN flow Elapse time: "
        //     << cuDFNsys::CPUSecond() - iStart_DFN << " seconds\n\n";
        TimeDFN_FLOW[i] = cuDFNsys::CPUSecond() - iStart_DFN;
    };

    chdir(curPath.c_str());
    cuDFNsys::HDF5API h5g;

    string FileNameH5 =
        "TimeElapsed_cuDFNsys_" + cuDFNsys::ToStringWithWidth(Num, 3) + ".h5";
    h5g.NewFile(FileNameH5);

    h5g.AddDataset<double>(FileNameH5, "N", "TimeDFN_GEN", TimeDFN_GEN.data(),
                           make_uint2(1, MCtime));
    h5g.AddDataset<double>(FileNameH5, "N", "TimeDFN_MESH", TimeDFN_MESH.data(),
                           make_uint2(1, MCtime));
    h5g.AddDataset<double>(FileNameH5, "N", "TimeDFN_FLOW", TimeDFN_FLOW.data(),
                           make_uint2(1, MCtime));

    return 0;
}