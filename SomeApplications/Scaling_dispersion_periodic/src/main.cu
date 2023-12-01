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

#include "cuDFNsys.cuh"
#include <H5Exception.h>
#include <fstream>
#include <iostream>
#include <limits.h>
#include <numeric>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>
string GetArgvWithoutComment(const string &A, const char &B);
int main(int argc, char *argv[])
{
    char cwd[PATH_MAX];
    if (getcwd(cwd, sizeof(cwd)) != NULL)
        printf("Current working dir: %s\n", cwd);
    else
        throw cuDFNsys::ExceptionsPause("getcwd() error");
    string curPath = cwd;
    //----------------------

    time_t t;
    time(&t);

    //----------------------
    char split_char = '#';

    uint DFN_No_init = atoi(GetArgvWithoutComment(argv[1], split_char).c_str());
    uint DFN_No_end = atoi(GetArgvWithoutComment(argv[2], split_char).c_str());
    uint InitNumFracs =
        atoi(GetArgvWithoutComment(argv[3], split_char).c_str());
    uint NumFracsIncreament =
        atoi(GetArgvWithoutComment(argv[4], split_char).c_str());
    uint NumTimeStepsPT_singleLoop =
        atoi(GetArgvWithoutComment(argv[5], split_char).c_str());
    uint NumLoops_PT = atoi(GetArgvWithoutComment(argv[6], split_char).c_str());

    double kappa_ = atof(GetArgvWithoutComment(argv[7], split_char).c_str());
    double DomainSize_x_ =
        atof(GetArgvWithoutComment(argv[8], split_char).c_str());
    double Beta_a = atof(GetArgvWithoutComment(argv[9], split_char).c_str());
    double Gamma_a = atof(GetArgvWithoutComment(argv[10], split_char).c_str());
    int ModeSizeDistri =
        atoi(GetArgvWithoutComment(argv[11], split_char).c_str());
    double4 SizeDistributionParameters_a =
        make_double4(atof(GetArgvWithoutComment(argv[12], split_char).c_str()),
                     atof(GetArgvWithoutComment(argv[13], split_char).c_str()),
                     atof(GetArgvWithoutComment(argv[14], split_char).c_str()),
                     atof(GetArgvWithoutComment(argv[15], split_char).c_str()));
    double minGridSize_a =
        atof(GetArgvWithoutComment(argv[16], split_char).c_str());
    double maxGridSize_a =
        atof(GetArgvWithoutComment(argv[17], split_char).c_str());

    double MuOverRhoG_a =
        atof(GetArgvWithoutComment(argv[18], split_char).c_str());
    double ConsTq_a = atof(GetArgvWithoutComment(argv[19], split_char).c_str());

    int NumRandomWalkers =
        atoi(GetArgvWithoutComment(argv[20], split_char).c_str());
    double Pe_a = atof(GetArgvWithoutComment(argv[21], split_char).c_str());
    double LengthScale_a =
        atof(GetArgvWithoutComment(argv[22], split_char).c_str());
    double FactorToDivideDeltaT =
        atof(GetArgvWithoutComment(argv[23], split_char).c_str());
    string injectionMode = GetArgvWithoutComment(argv[24], split_char).c_str();
    double customedInjectionPlane =
        atof(GetArgvWithoutComment(argv[25], split_char).c_str());
    bool IfUseFluxWeightedOrEqualProbableMixingIntersection_a =
        (atoi(GetArgvWithoutComment(argv[26], split_char).c_str()) == 0 ? false
                                                                        : true);
    bool If_ReRun_ =
        (atoi(GetArgvWithoutComment(argv[27], split_char).c_str()) == 0 ? false
                                                                        : true);
    int ReRunMaxTime =
        atoi(GetArgvWithoutComment(argv[28], split_char).c_str());
    string LogFile = GetArgvWithoutComment(argv[29], split_char).c_str();

    //------------------------------
    std::ofstream oss("./InputInfo.info", ios::out);
    oss << "DFN_No_init: " << DFN_No_init << endl;
    oss << "DFN_No_end: " << DFN_No_end << endl;
    oss << "InitNumFracs: " << InitNumFracs << endl;
    oss << "NumFracsIncreament: " << NumFracsIncreament << endl;
    oss << "NumTimeStepsPT_singleLoop: " << NumTimeStepsPT_singleLoop << endl;
    oss << "NumLoops_PT: " << NumLoops_PT << endl;
    oss << "kappa_: " << kappa_ << endl;
    oss << "DomainSize_x_: " << DomainSize_x_ << endl;
    oss << "Beta_a: " << Beta_a << endl;
    oss << "Gamma_a: " << Gamma_a << endl;
    oss << "ModeSizeDistri: " << ModeSizeDistri << endl;
    oss << "SizeDistributionParameters_a: " << SizeDistributionParameters_a.x
        << ", " << SizeDistributionParameters_a.y << ", "
        << SizeDistributionParameters_a.z << ", "
        << SizeDistributionParameters_a.w << endl;
    oss << "minGridSize_a: " << minGridSize_a << endl;
    oss << "maxGridSize_a: " << maxGridSize_a << endl;
    oss << "MuOverRhoG_a: " << MuOverRhoG_a << endl;
    oss << "ConsTq_a: " << ConsTq_a << endl;
    oss << "NumRandomWalkers: " << NumRandomWalkers << endl;
    oss << "Pe_a: " << Pe_a << endl;
    oss << "LengthScale_a: " << LengthScale_a << endl;
    oss << "FactorToDivideDeltaT: " << FactorToDivideDeltaT << endl;
    oss << "injectionMode: " << injectionMode << endl;
    oss << "customedInjectionPlane: " << customedInjectionPlane << endl;
    oss << "IfUseFluxWeightedOrEqualProbableMixingIntersection_a: "
        << IfUseFluxWeightedOrEqualProbableMixingIntersection_a << endl;
    oss << "If_ReRun_: " << If_ReRun_ << endl;
    oss << "ReRunMaxTime: " << ReRunMaxTime << endl;
    oss << "LogFile: " << LogFile << endl;
    oss.close();
    //----------------------
    if (DFN_No_init <= 0)
    {
        cout << "`DFN_No_init` cannot be <= 0\n";
        return 0;
    }
    if (DFN_No_end < DFN_No_init)
    {
        cout << "`DFN_No_end` cannot be < `DFN_No_init`\n";
        return 0;
    }

    int dev = 0;
    GPUErrCheck(cudaSetDevice(dev));
    cuDFNsys::Warmup<<<256 / 256 + 1, 256 /*  1, 2*/>>>();
    cudaDeviceSynchronize();

    //----------------------
    cuDFNsys::HDF5API h5g;
    for (uint i = DFN_No_init; i <= DFN_No_end; ++i)
    {
        string path2 = "DFN_" + cuDFNsys::ToStringWithWidth(i, 3);
        string command1 = "mkdir -p " + path2;
        system(command1.c_str());

        string command2 = curPath + "/" + path2;
        chdir(command2.c_str());

        int NumFracs = InitNumFracs + (i - 1) * NumFracsIncreament;

        if (i == DFN_No_init)
            system("echo \" \" > ../recordTime_and_Error.log");

        string Start_i = "loop " + std::to_string(i) +
                         ", NF = " + std::to_string(NumFracs) + " started\n";
        system("echo $(date +%d-%m-%y---%T) >> ../recordTime_and_Error.log");
        Start_i = "echo \"" + Start_i + "\" >> ../recordTime_and_Error.log";
        system(Start_i.c_str());

        uint countReRunTime = 0; // ReRunMaxTime

        cuDFNsys::DFN<double> my_dfn;
        cuDFNsys::MeshDFN<double> meshGen;
        cuDFNsys::FlowDFN<double> flowDFN;
        cuDFNsys::PTDFN<double> particleTracking;

        for (int j = 0; j < NumLoops_PT; ++j)
        {
            try
            {
                string GGF = "echo \" \" > ../" + LogFile;

                system(GGF.c_str());

                cout << "NumFracs: " << NumFracs << ", j = " << j << endl;

                std::ifstream File_d1("./Class_DFN.h5"),
                    File_d2("./Class_MESH.h5"), File_d3("./Class_FLOW.h5");

                bool DFSW =
                    ((File_d1.good() && File_d2.good() && File_d3.good())
                         ? true
                         : false);
                File_d1.close(), File_d2.close(), File_d3.close();
                if (!DFSW)
                {
                    for (int k = 0; k < 100000; ++k)
                    {
                        my_dfn.RandomSeed = (unsigned long)(t + k * 1e5);
                        my_dfn.NumFractures = {NumFracs};
                        my_dfn.Kappa = {kappa_};
                        my_dfn.MeanOrientationOfFisherDistribution = {
                            make_double3(0., 0., 1.)};
                        my_dfn.DomainSizeX = DomainSize_x_;
                        my_dfn.DomainDimensionRatio = make_double3(1., 1., 1.);
                        my_dfn.Beta = {Beta_a};
                        my_dfn.Gamma = {Gamma_a};
                        my_dfn.ModeOfSizeDistribution = {ModeSizeDistri};
                        my_dfn.SizeDistributionParameters = {
                            SizeDistributionParameters_a};
                        my_dfn.PercoDir = 2;
                        my_dfn.FractureGeneration();
                        my_dfn.IdentifyIntersectionsClusters(true);
                        if (my_dfn.PercolationCluster.size() > 0)
                            break;
                        if (k == 99999)
                        {
                            cout << "Did not find percolation cluster\n";
                            return 0;
                        }
                    }

                    if (my_dfn.PercolationCluster.size() > 0)
                    {
                        my_dfn.IdentifyIntersectionsClusters(false);
                        std::vector<double> Data1(14);
                        cuDFNsys::GetStatistics<double>(
                            my_dfn.FracturesHost, my_dfn.IntersectionMap,
                            my_dfn.ListClusters, my_dfn.PercolationCluster,
                            my_dfn.DomainSizeX, Data1[0], Data1[1], Data1[2],
                            Data1[3], Data1[4], Data1[5], Data1[6], Data1[7],
                            Data1[8], Data1[9], Data1[10], Data1[11], Data1[12],
                            Data1[13]);

                        my_dfn.SpatialPeriodicity();
                        my_dfn.IdentifyIntersectionsClusters(true);
                        my_dfn.Visualization("DFN_VISUAL_I", "DFN_VISUAL_I",
                                             "DFN_VISUAL_I", true, true, true,
                                             true);

                        meshGen.MinElementSize = minGridSize_a;
                        meshGen.MaxElementSize = maxGridSize_a;
                        meshGen.MeshGeneration(my_dfn);
                        meshGen.Visualization(my_dfn, "DFN_MESH_VISUAL",
                                              "DFN_MESH_VISUAL",
                                              "DFN_MESH_VISUAL", true, true);

                        flowDFN.IfPeriodic = true;
                        flowDFN.MuOverRhoG = MuOverRhoG_a;
                        flowDFN.ConsTq = ConsTq_a;

                        flowDFN.FlowSimulation(my_dfn, meshGen);
                        flowDFN.Visualization(
                            my_dfn, meshGen, "DFN_FLOW_VISUAL",
                            "DFN_FLOW_VISUAL", "DFN_FLOW_VISUAL");
                        flowDFN.StoreInH5("Class_FLOW");
                        meshGen.StoreInH5("Class_MESH");
                        my_dfn.StoreInH5("Class_DFN");

                        h5g.NewFile("FlowProperties.h5");
                        vector<string> datasetname = {
                            "P33_total",
                            "P33_connected_z",
                            "Ratio_of_P33_z",
                            "P33_largest_cluster",
                            "P32_total",
                            "P32_connected_z",
                            "Ratio_of_P32_z",
                            "P32_largest_cluster",
                            "P30",
                            "P30_connected_z",
                            "Ratio_of_P30_z",
                            "P30_largest_cluster",
                            "Percolation_probability_z",
                            "n_I",
                            "Permeability_z",
                            "Q_error_z"};
                        vector<double *> data_input = {
                            &Data1[0],
                            &Data1[1],
                            &Data1[2],
                            &Data1[3],
                            &Data1[4],
                            &Data1[5],
                            &Data1[6],
                            &Data1[7],
                            &Data1[8],
                            &Data1[9],
                            &Data1[10],
                            &Data1[11],
                            &Data1[12],
                            &Data1[13],
                            &flowDFN.FlowData.Permeability,
                            &flowDFN.FlowData.QError};
                        vector<uint2> dim_ss(data_input.size(),
                                             make_uint2(1, 1));
                        h5g.AddDatasetsWithOneGroup<double>(
                            "FlowProperties.h5", "ConnectivityPermeability",
                            datasetname, data_input, dim_ss);
                        h5g.AddDataset<double>("FlowProperties.h5", "N",
                                               "MeanV", &flowDFN.MeanVelocity,
                                               make_uint2(1, 1));
                        h5g.AddDataset<double>("FlowProperties.h5", "N", "MaxV",
                                               &flowDFN.MaxVelocity,
                                               make_uint2(1, 1));
                    }
                }
                else
                {
                    my_dfn.LoadClassFromH5("Class_DFN");
                    my_dfn.IdentifyIntersectionsClusters(true);
                    meshGen.LoadClassFromH5("Class_MESH");
                    flowDFN.LoadClassFromH5("Class_FLOW");

                    std::ifstream File_d6("../scripts/DelDataSet.py");
                    if (!File_d6.good())
                    {
                        File_d6.close();
                        cout << "`scripts/DelDataSet.py` is not existing\n";
                        return 0;
                    }
                    system("python3 ../scripts/DelDataSet.py");
                }
                particleTracking.IfPeriodic = true;
                particleTracking.NumTimeSteps = NumTimeStepsPT_singleLoop;
                particleTracking.NumParticles = NumRandomWalkers;
                particleTracking.MolecularDiffusion =
                    LengthScale_a * flowDFN.FlowData.MeanVelocity / Pe_a;

                //-----------------------------------advction time to cross a grid
                double DeltaT_Advection =
                    pow(meshGen.MeshData.MeanGridSize, 0.5) /
                    flowDFN.FlowData.MaxVelocity;
                DeltaT_Advection /= FactorToDivideDeltaT;
                //-----------------------------------diffusion time to cross a grid
                double DeltaT_Diffusion =
                    pow(meshGen.MeshData.MeanGridSize, 0.5) /
                    (pow(2 * particleTracking.MolecularDiffusion, 0.5) * 4);
                DeltaT_Diffusion =
                    DeltaT_Diffusion * DeltaT_Diffusion / FactorToDivideDeltaT;
                cout << "DeltaT_Advection: " << DeltaT_Advection
                     << ", DeltaT_Diffusion: " << DeltaT_Diffusion << endl;
                if (DeltaT_Advection < DeltaT_Diffusion)
                {
                    cout << "\n** ADVEction Time is smaller **\n\n";
                    particleTracking.DeltaT = DeltaT_Advection;
                }
                else
                {
                    cout << "\n** DIFFusion Time is smaller **\n\n";
                    particleTracking.DeltaT = DeltaT_Diffusion;
                }
                particleTracking.InjectionMethod = injectionMode;
                particleTracking.OutputAllPTInformationOrFPTCurve = false;
                particleTracking.SpacingOfControlPlanes = 3000000000;
                particleTracking.IfOutputVarianceOfDisplacementsEachStep = true;
                particleTracking.IfInjectAtCustomedPlane = true;
                particleTracking.CustomedPlaneInjection =
                    customedInjectionPlane;
                particleTracking
                    .IfUseFluxWeightedOrEqualProbableMixingIntersection =
                    IfUseFluxWeightedOrEqualProbableMixingIntersection_a;
                particleTracking.ParticleTracking(my_dfn, meshGen, flowDFN);
            }
            catch (cuDFNsys::ExceptionsIgnore &e)
            {
                cout << "cuDFNsys::ExceptionsIgnore\n";
                cout << e.what() << endl;
                cout << path2 << endl;
                system("rm -rf *.h5 ParticlePositionResult");
                j = 0;

                system("echo $(date +%d-%m-%y---%T) >> "
                       "../recordTime_and_Error.log");
                string ERDS = "echo \"\t" + e.what() +
                              "\n\n\" >> ../recordTime_and_Error.log";
                system(ERDS.c_str());
            }
            catch (cuDFNsys::ExceptionsPause &e)
            {
                cout << "cuDFNsys::ExceptionsPause\n";
                cout << e.what() << endl;
                cout << path2 << endl;
                if (!If_ReRun_)
                {
                    system("rm -rf *.h5 ParticlePositionResult");
                    j = 0;
                }
                else
                { // re-run
                    string vccy =
                        "result=$(grep \"Loop times is too large!\" ../" +
                        LogFile + ")";
                    int status = system(vccy.c_str());

                    if (status == 0)
                    {
                        countReRunTime++;

                        if (countReRunTime <= ReRunMaxTime)
                        {
                            system("echo $(date +%d-%m-%y---%T) >> "
                                   "../recordTime_and_Error.log");
                            string ERDS = "echo \"ReRunTheSameDFN_Transport\" "
                                          ">> ../recordTime_and_Error.log";
                            system(ERDS.c_str());
                        }
                        else
                        {
                            string ERDS =
                                "echo \"ReRunTimes reached the maximum\" >> "
                                "../recordTime_and_Error.log";
                            system(ERDS.c_str());

                            system("rm -rf *.h5 ParticlePositionResult");
                            j = 0;
                        }
                        //exit(0);
                    }
                    else
                    {
                        system("rm -rf *.h5 ParticlePositionResult");
                        j = 0;
                        //exit(0);
                    }
                }

                system("echo $(date +%d-%m-%y---%T) >> "
                       "../recordTime_and_Error.log");
                string ERDS = "echo \"\t" + e.what() +
                              "\n\n\" >> ../recordTime_and_Error.log";
                system(ERDS.c_str());
            }
            catch (H5::Exception &e)
            {
                cout << "H5::Exception\n";
                //e.printError();
                cout << path2 << endl;
                system("rm -rf *.h5 ParticlePositionResult");
                j = 0;

                system("echo $(date +%d-%m-%y---%T) >> "
                       "../recordTime_and_Error.log");
                system("echo \"H5::Exception\n\n\" >> "
                       "../recordTime_and_Error.log");
            }
            catch (...)
            {
                cout << "Unknown exceptions!\n";
                cout << path2 << endl;
                system("rm -rf *.h5 ParticlePositionResult");
                j = 0;

                system("echo $(date +%d-%m-%y---%T) >> "
                       "../recordTime_and_Error.log");
                system("echo \"Unknown exceptions\n\n\" >> "
                       "../recordTime_and_Error.log");
            }
        }
        chdir(curPath.c_str());
    }

    return 0;
};

string GetArgvWithoutComment(const string &A, const char &B)
{
    std::stringstream test(A);
    std::string segment;
    std::vector<std::string> seglist;

    while (std::getline(test, segment, B))
    {
        seglist.push_back(segment);
    }

    return seglist[0];
};