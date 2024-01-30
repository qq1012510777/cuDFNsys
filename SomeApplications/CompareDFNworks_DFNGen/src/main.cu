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

    char cwd[PATH_MAX];
    if (getcwd(cwd, sizeof(cwd)) != NULL)
        printf("Current working dir: %s\n", cwd);
    else
        throw cuDFNsys::ExceptionsPause("getcwd() error");
    string curPath = cwd;

    int NumFractures_init = 1500;
    int FracIncreament = 200;

    int Num_MCTimes = 30;
    int Num_FracIncrements = 10;

    std::vector<std::vector<double>> TimeDFN_GEN(Num_FracIncrements);

    for (int i = 0; i < Num_FracIncrements; ++i)
    {
        TimeDFN_GEN[i].resize(Num_MCTimes);

        for (int j = 0; j < Num_MCTimes; ++j)
        {
            chdir(curPath.c_str());
            string path2 = "DFN_NumFrac_" +
                           cuDFNsys::ToStringWithWidth(
                               NumFractures_init + i * FracIncreament, 4) +
                           "_MC_NO_" + cuDFNsys::ToStringWithWidth(j + 1, 3);
            string command1 = "mkdir -p " + path2;
            system(command1.c_str());

            string command2 = curPath + "/" + path2;
            chdir(command2.c_str());

            cuDFNsys::DFN<double> my_dfn;

            my_dfn.NumFractures = {NumFractures_init + i * FracIncreament};
            my_dfn.Kappa = {0};
            my_dfn.MeanOrientationOfFisherDistribution = {
                make_double3(0, 0, 1)};
            my_dfn.DomainSizeX = 100;
            my_dfn.DomainDimensionRatio = make_double3(1, 1, 1);
            my_dfn.Beta = {0.3};
            my_dfn.Gamma = {1e-10};
            my_dfn.ModeOfSizeDistribution = {3};
            my_dfn.SizeDistributionParameters = {make_double4(7.5, 0, 0, 0)};
            my_dfn.PercoDir = 2;
            my_dfn.RandomSeed = (unsigned long)t + j;

            double iStart_DFN = cuDFNsys::CPUSecond();
            my_dfn.FractureGeneration();
            my_dfn.IdentifyIntersectionsClusters(true);
            TimeDFN_GEN[i][j] = cuDFNsys::CPUSecond() - iStart_DFN;

            my_dfn.StoreInH5("Class_DFN");
            my_dfn.Visualization("DFN_VISUAL", "DFN_VISUAL", "DFN_VISUAL", true,
                                 true, true, true);
        };
    };

    chdir(curPath.c_str());
    cuDFNsys::HDF5API h5g;
    h5g.NewFile("TimeElapsed_cuDFNsys_DFNGen.h5");

    for (int i = 0; i < Num_FracIncrements; ++i)
        h5g.AddDataset<double>(
            "TimeElapsed_cuDFNsys_DFNGen.h5", "N",
            "TimeDFN_GEN_NumFrac_" +
                cuDFNsys::ToStringWithWidth(
                    NumFractures_init + i * FracIncreament, 4),
            TimeDFN_GEN[i].data(), make_uint2(1, Num_MCTimes));

    return 0;
}