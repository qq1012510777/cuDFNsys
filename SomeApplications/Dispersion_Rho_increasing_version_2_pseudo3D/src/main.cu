#include "cuDFNsys.cuh"
#include <algorithm>
#include <chrono>
#include <cstdlib> // For std::system
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits.h>
#include <memory>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <unistd.h>
#include <vector>
using namespace std;
namespace fs = std::filesystem;

bool RunCMD_with_RealTimeCheck(const string &cmd, const string &logFile,
                               const bool &IfGenLogFile = false);
void CreateOrEmptyFile(const std::string &filename);
void AddLineToFile(const std::string &filename, const std::string &line);
bool IfAFileExist(const string &FilePath);
string DoubleNumberToScientificNotationString(const double &number);
uint GetH5DatasetSize(const string &nameH5, const string &nameDataset);

int main(int argc, char *argv[])
{
    time_t t;
    time(&t);

    string ExeuctablePath = argv[1];
    // DFN generation parameters
    double DomainSizeX = atof(argv[2]);
    double3 DomainDimensionRatio =
        make_double3(atof(argv[3]), atof(argv[4]), atof(argv[5]));
    int PercoDir = atoi(argv[6]);
    double Kappa = atof(argv[7]);

    int FracNumInit = atoi(argv[8]);
    int FracNumIncre = atoi(argv[9]);
    int NumFracIncre = atoi(argv[10]);

    int FracSizeDistriMode = atoi(argv[11]);
    double4 FracSizeDistriPara = make_double4(atof(argv[12]), atof(argv[13]),
                                              atof(argv[14]), atof(argv[15]));

    double Beta = atof(argv[16]);
    double Gamma = atof(argv[17]);
    // DFN mesh parameters
    double MeshMinimumGridSize = atof(argv[18]);
    double MeshMaximumGridSize = atof(argv[19]);

    // DFN PT parameters
    int NumberParticles = atoi(argv[20]);
    int NumSteps = atoi(argv[21]);
    string InjectionMethod_initialcondition = argv[22];
    int If_OutputAllPTInformationOrFPTCurve = atoi(argv[23]);
    double SpacingOfControlPlanes = atof(argv[24]);
    int IfOutputVarianceOfDisplacementsEachStep = atoi(argv[25]);
    int IfInjectAtCustomedPlane = atoi(argv[26]);
    double CustomedPlaneInjection = atof(argv[27]);
    int IfUseFluxWeightedOrEqualProbableMixingIntersection = atoi(argv[28]);
    int TimeIntervalOutPTInformation = atoi(argv[29]);
    int IfOutputAllParticlesAccumulativeDisplacements = atoi(argv[30]);

    ///
    double Pe = atof(argv[31]);
    double LengthScalePe = atof(argv[32]);
    double FactorReduceCharacteristicTimescale = atof(argv[33]);

    /// termination condition for PT
    double Percentage_Particle_escape_from_inlet = atof(argv[34]);
    double Percentage_Particle_arrived_outlet = atof(argv[35]);
    uint ErrorCountAllowable = atoi(argv[36]);
    uint NumStepsThreshold = atoi(argv[37]);

    cout << "NumberParticles: " << NumberParticles << endl;
    cout << "NumSteps: " << NumSteps << endl;
    cout << "InjectionMethod_initialcondition: "
         << InjectionMethod_initialcondition << endl;
    cout << "If_OutputAllPTInformationOrFPTCurve: "
         << If_OutputAllPTInformationOrFPTCurve << endl;
    cout << "SpacingOfControlPlanes: " << SpacingOfControlPlanes << endl;
    cout << "IfOutputVarianceOfDisplacementsEachStep: "
         << IfOutputVarianceOfDisplacementsEachStep << endl;
    cout << "IfInjectAtCustomedPlane: " << IfInjectAtCustomedPlane << endl;
    cout << "CustomedPlaneInjection: " << CustomedPlaneInjection << endl;
    cout << "IfUseFluxWeightedOrEqualProbableMixingIntersection: "
         << IfUseFluxWeightedOrEqualProbableMixingIntersection << endl;
    cout << "TimeIntervalOutPTInformation: " << TimeIntervalOutPTInformation
         << endl;
    cout << "IfOutputAllParticlesAccumulativeDisplacements: "
         << IfOutputAllParticlesAccumulativeDisplacements << endl;
    cout << "Pe: " << Pe << endl;
    cout << "LengthScalePe: " << LengthScalePe << endl;
    cout << "FactorReduceCharacteristicTimescale: "
         << FactorReduceCharacteristicTimescale << endl;

    cout << "If_OutputAllPTInformationOrFPTCurve: "
         << If_OutputAllPTInformationOrFPTCurve << endl;
    cout << "Percentage_Particle_escape_from_inlet: "
         << Percentage_Particle_escape_from_inlet << endl;
    cout << "Percentage_Particle_arrived_outlet: "
         << Percentage_Particle_arrived_outlet << endl;
    cout << "ErrorCountAllowable: " << ErrorCountAllowable << endl;
    cout << "NumStepsThreshold: " << NumStepsThreshold << endl;
    //exit(0);
    cuDFNsys::HDF5API h5g;

    for (int i = 0; i <= NumFracIncre; ++i)
    {

        string DFNFileName = "DFN_" + cuDFNsys::ToStringWithWidth(i + 1, 3);

        int result_system = system(("mkdir ./" + DFNFileName + " -p").c_str());

        string LogFile =
            "log_DFN_" + cuDFNsys::ToStringWithWidth(i + 1, 3) + ".txt";

        CreateOrEmptyFile("./" + DFNFileName + "/" + LogFile);
        try
        {
            //-----------------DFN_Gen-------------
            if (!IfAFileExist("./" + DFNFileName + "/Class_DFN.h5"))
            {
                string DFN_Gen_csv_name = "StochasticDFN";

                CreateOrEmptyFile(DFNFileName + "/" + DFN_Gen_csv_name +
                                  ".csv");

                AddLineToFile("./" + DFNFileName + "/" + DFN_Gen_csv_name +
                                  ".csv",
                              "IfStochastic,1,\n");
                AddLineToFile(
                    "./" + DFNFileName + "/" + DFN_Gen_csv_name + ".csv",
                    "DomainSizeX," + std::to_string(DomainSizeX) + ",\n");
                AddLineToFile(
                    "./" + DFNFileName + "/" + DFN_Gen_csv_name + ".csv",
                    "DomainDimensionRatio," +
                        std::to_string(DomainDimensionRatio.x) + "," +
                        std::to_string(DomainDimensionRatio.y) + "," +
                        std::to_string(DomainDimensionRatio.z) + ",\n");
                AddLineToFile("./" + DFNFileName + "/" + DFN_Gen_csv_name +
                                  ".csv",
                              "Percolation_direction," +
                                  std::to_string(PercoDir) + ",\n");
                AddLineToFile("./" + DFNFileName + "/" + DFN_Gen_csv_name +
                                  ".csv",
                              "NumFractureGroups,1,\n");
                AddLineToFile(
                    "./" + DFNFileName + "/" + DFN_Gen_csv_name + ".csv",
                    "NumFractureEachGroup," +
                        std::to_string(FracNumInit + i * FracNumIncre) + ",\n");
                AddLineToFile("./" + DFNFileName + "/" + DFN_Gen_csv_name +
                                  ".csv",
                              "KappaValues," + std::to_string(Kappa) + ",\n");
                AddLineToFile("./" + DFNFileName + "/" + DFN_Gen_csv_name +
                                  ".csv",
                              "MeanOrientationOfFisherDistribution,0,0,1,\n");
                AddLineToFile("./" + DFNFileName + "/" + DFN_Gen_csv_name +
                                  ".csv",
                              "ModeOfSizeDistribution," +
                                  std::to_string(FracSizeDistriMode) + ",\n");
                AddLineToFile("./" + DFNFileName + "/" + DFN_Gen_csv_name +
                                  ".csv",
                              "SizeDistributionParameters," +
                                  std::to_string(FracSizeDistriPara.x) + "," +
                                  std::to_string(FracSizeDistriPara.y) + "," +
                                  std::to_string(FracSizeDistriPara.z) + "," +
                                  std::to_string(FracSizeDistriPara.w) + ",\n");
                AddLineToFile("./" + DFNFileName + "/" + DFN_Gen_csv_name +
                                  ".csv",
                              "Beta," + std::to_string(Beta) + ",\n");
                AddLineToFile(
                    "./" + DFNFileName + "/" + DFN_Gen_csv_name + ".csv",
                    "Gamma," + DoubleNumberToScientificNotationString(Gamma) +
                        ",\n");
                string DFN_gen_run_command = "cd ./" + DFNFileName + " && " +
                                             ExeuctablePath + "/DFN_Gen " +
                                             "./StochasticDFN";

                for (int j = 0; j < 10000; ++j)
                {
                    bool DFN_gen_success = RunCMD_with_RealTimeCheck(
                        DFN_gen_run_command, "./" + DFNFileName + "/" + LogFile,
                        true);
                    int PercolationClustersSize =
                        GetH5DatasetSize("./" + DFNFileName + "/Class_DFN.h5",
                                         "PercolationClusters");

                    if (!DFN_gen_success || PercolationClustersSize == 0)
                    {
                        if (!DFN_gen_success)
                            cout << DFNFileName << ": DFN_Gen failed\n";
                        result_system =
                            system(("rm -rf ./" + DFNFileName + "/Class_DFN.h5")
                                       .c_str());
                        continue;
                    }

                    if (PercolationClustersSize > 0)
                    {
                        system(("cp ./" + DFNFileName + "/Class_DFN.h5  ./" +
                                DFNFileName + "/Class_DFN_original.h5")
                                   .c_str());
                        break;
                    }
                }
            }

            //-----------------DFN_Mesh------------
            if (IfAFileExist("./" + DFNFileName + "/Class_DFN.h5") &&
                !IfAFileExist("./" + DFNFileName + "/Class_MESH.h5"))
            {
                string DFN_Mesh_csv_name = "MeshPara";

                CreateOrEmptyFile(DFNFileName + "/" + DFN_Mesh_csv_name +
                                  ".csv");

                AddLineToFile("./" + DFNFileName + "/" + DFN_Mesh_csv_name +
                                  ".csv",
                              "ExpectedMinimimGridSize," +
                                  std::to_string(MeshMinimumGridSize) + ",\n");
                AddLineToFile("./" + DFNFileName + "/" + DFN_Mesh_csv_name +
                                  ".csv",
                              "ExpectedMaximimGridSize," +
                                  std::to_string(MeshMaximumGridSize) + ",\n");
                string DFN_mesh_run_command = "cd ./" + DFNFileName + " && " +
                                              ExeuctablePath + "/DFN_Mesh " +
                                              "./MeshPara";
                bool DFN_mesh_success = RunCMD_with_RealTimeCheck(
                    DFN_mesh_run_command, "./" + DFNFileName + "/" + LogFile);
                if (!DFN_mesh_success)
                {
                    if (!DFN_mesh_success)
                        cout << DFNFileName
                             << ": DFN_Mesh failed, regenerate DFN\n";
                    result_system =
                        system(("rm -rf ./" + DFNFileName + "/Class_DFN.h5 ./" +
                                DFNFileName + "/Class_MESH.h5")
                                   .c_str());
                    i--;
                    continue;
                }
            }

            //-----------------DFN_Flow---------------
            if (IfAFileExist("./" + DFNFileName + "/Class_DFN.h5") &&
                IfAFileExist("./" + DFNFileName + "/Class_MESH.h5") &&
                !IfAFileExist("./" + DFNFileName + "/Class_FLOW.h5"))
            {
                string DFN_flow_csv_name = "FlowPara";

                CreateOrEmptyFile(DFNFileName + "/" + DFN_flow_csv_name +
                                  ".csv");

                AddLineToFile("./" + DFNFileName + "/" + DFN_flow_csv_name +
                                  ".csv",
                              "MuOverRhoG,1,\n");
                AddLineToFile(
                    "./" + DFNFileName + "/" + DFN_flow_csv_name + ".csv",
                    "InletHead," + std::to_string(DomainSizeX) + ",\n");
                AddLineToFile("./" + DFNFileName + "/" + DFN_flow_csv_name +
                                  ".csv",
                              "OutletHead,0,\n");
                string DFN_flow_run_command = "cd ./" + DFNFileName + " && " +
                                              ExeuctablePath + "/DFN_Flow " +
                                              "./FlowPara";
                bool DFN_flow_success = RunCMD_with_RealTimeCheck(
                    DFN_flow_run_command, "./" + DFNFileName + "/" + LogFile);
                if (!DFN_flow_success)
                {
                    if (!DFN_flow_success)
                        cout << DFNFileName
                             << ": DFN_Flow failed, regenerate DFN\n";
                    result_system =
                        system(("rm -rf ./" + DFNFileName + "/Class_DFN.h5 ./" +
                                DFNFileName + "/Class_MESH.h5 ./" +
                                DFNFileName + "/Class_FLOW.h5")
                                   .c_str());
                    i--;
                    continue;
                }
            }

            //-----------------DFN_PT---------------
            if (IfAFileExist("./" + DFNFileName + "/Class_DFN.h5") &&
                IfAFileExist("./" + DFNFileName + "/Class_MESH.h5") &&
                IfAFileExist("./" + DFNFileName + "/Class_FLOW.h5"))
            {
                std::vector<double> tmp;

                // check if PT is normal
                if (IfAFileExist("./" + DFNFileName +
                                 "/ParticlePositionResult/DispersionInfo.h5"))
                {
                    result_system = system(
                        ("rm -rf ./" + DFNFileName + "/PTFinished").c_str());
                    cout << DFNFileName
                         << ": checking particles escaping from inlet and "
                            "outlet\n";
                    tmp = h5g.ReadDataset<double>(
                        "./" + DFNFileName +
                            "/ParticlePositionResult/DispersionInfo.h5",
                        "N", "NumParticles");
                    int NumParticlesTotal = tmp[0];

                    tmp = h5g.ReadDataset<double>(
                        "./" + DFNFileName +
                            "/ParticlePositionResult/DispersionInfo.h5",
                        "N", "NumParticlesLeftFromInlet");
                    int NumParticlesLeftFromInlet = tmp[0];

                    /// if too many particles escaped from inlet
                    if (NumParticlesLeftFromInlet * 1.0 /
                            (NumParticlesTotal * 1.0) >
                        Percentage_Particle_escape_from_inlet)
                    {
                        cout << DFNFileName
                             << ": DFN_PT failed. Too many particles escaped "
                                "from "
                                "the "
                                "inlet. Regenerate DFN\n";
                        result_system = system(("rm -rf ./" + DFNFileName +
                                                "/*.h5 " + DFNFileName +
                                                "/ParticlePositionResult/* ./" +
                                                DFNFileName + "/PTFinished")
                                                   .c_str());
                        i--;
                        continue;
                    }

                    // if almost all particles arrived
                    std::vector<uint> FPT = h5g.ReadDataset<uint>(
                        "./" + DFNFileName +
                            "/ParticlePositionResult/"
                            "ParticlePosition_FPTControlPlanes.h5",
                        "N",
                        "ControlPlane_" +
                            cuDFNsys::ToStringWithWidth(DomainSizeX, 10) +
                            "_m");
                    int zeroCount = std::count(FPT.begin(), FPT.end(), 0);
                    if ((NumParticlesTotal - zeroCount) * 1.0 /
                            (NumParticlesTotal * 1.0) >
                        Percentage_Particle_arrived_outlet)
                    {
                        cout << DFNFileName << ": DFN_PT finished, "
                             << NumParticlesTotal - zeroCount
                             << " particles arrived\n";
                        CreateOrEmptyFile(DFNFileName + "/PTFinished");
                        continue;
                    }

                    // if too many steps have been run
                    std::vector<uint> tmp_uinty = h5g.ReadDataset<uint>(
                        "./" + DFNFileName +
                            "/ParticlePositionResult/DispersionInfo.h5",
                        "N", "NumOfSteps");
                    cout << DFNFileName << ": " << NumParticlesTotal - zeroCount
                         << "/" << NumParticlesTotal
                         << " arrived, NumOfSteps: " << tmp_uinty[0] << endl;
                    if (tmp_uinty[0] > NumStepsThreshold)
                    {
                        cout << DFNFileName
                             << ": it runs too many steps:  " << tmp_uinty[0]
                             << "; " << zeroCount << "/" << NumParticlesTotal
                             << " has arrived\n";

                        //CreateOrEmptyFile(DFNFileName + "/PTFinished");
                        result_system = system(("rm -rf ./" + DFNFileName +
                                                "/*.h5 " + DFNFileName +
                                                "/ParticlePositionResult/* ./" +
                                                DFNFileName + "/PTFinished")
                                                   .c_str());
                        cout << "*********************************\n";
                        cout << "The current DFN is deleted\n";
                        cout << "*********************************\n";
                        i--;
                        continue;
                    }
                } //

                // if too many error happens during PT

                if (IfAFileExist("./" + DFNFileName + "/PT_ErrorCount.h5"))
                {
                    cout << DFNFileName
                         << ": checking the number of errors of PT\n";
                    std::vector<uint> tmp_count = h5g.ReadDataset<uint>(
                        "./" + DFNFileName + "/PT_ErrorCount.h5", "N",
                        "ErrorCount");
                    if (tmp_count[0] > ErrorCountAllowable)
                    {
                        cout << DFNFileName
                             << ": DFN_PT failed. Too many errors happened "
                                "during "
                                "PT. Regenerate DFN\n";
                        result_system = system(("rm -rf ./" + DFNFileName +
                                                "/*.h5 " + DFNFileName +
                                                "/ParticlePositionResult/* ./" +
                                                DFNFileName + "/PTFinished")
                                                   .c_str());
                        i--;
                        continue;
                    }
                }

                //     write PT csv
                {

                    tmp = h5g.ReadDataset<double>("./" + DFNFileName +
                                                      "/Class_FLOW.h5",
                                                  "N", "MeanVelocity");
                    double meanV = tmp[0];

                    tmp = h5g.ReadDataset<double>("./" + DFNFileName +
                                                      "/Class_FLOW.h5",
                                                  "N", "MaxVelocity");
                    double maxV = tmp[0];

                    tmp = h5g.ReadDataset<double>("./" + DFNFileName +
                                                      "/Class_MESH.h5",
                                                  "N", "MeanGridSize");
                    double meanGridSize = tmp[0];

                    double AdvectionTimeScale = pow(meanGridSize, 0.5) / maxV;
                    // FactorReduceCharacteristicTimescale
                    double MolecularDiffusionCoefficient =
                        meanV * LengthScalePe / Pe;

                    double DeltaT = 0;
                    cout << "MolecularDiffusionCoefficient = "
                         << MolecularDiffusionCoefficient << endl;
                    if (Pe == 0 || LengthScalePe == 0 ||
                        MolecularDiffusionCoefficient == 0)
                    {
                        MolecularDiffusionCoefficient = 0;

                        DeltaT = AdvectionTimeScale /
                                 FactorReduceCharacteristicTimescale;

                        cout << DFNFileName << ": advection only, DeltaT = "
                             << AdvectionTimeScale << " / "
                             << FactorReduceCharacteristicTimescale << " = "
                             << DeltaT << endl;
                    }
                    else
                    {
                        double DiffusionTimeScale =
                            meanGridSize / (2 * MolecularDiffusionCoefficient);
                        if (DiffusionTimeScale < AdvectionTimeScale)
                        {
                            DeltaT = DiffusionTimeScale /
                                     FactorReduceCharacteristicTimescale;
                            cout << DFNFileName
                                 << ": diffusion dominates, DeltaT = "
                                 << DiffusionTimeScale << " / "
                                 << FactorReduceCharacteristicTimescale << " = "
                                 << DeltaT << endl;
                        }
                        else
                        {
                            DeltaT = AdvectionTimeScale /
                                     FactorReduceCharacteristicTimescale;
                            cout << DFNFileName
                                 << ": advection dominates, DeltaT = "
                                 << AdvectionTimeScale << " / "
                                 << FactorReduceCharacteristicTimescale << " = "
                                 << DeltaT << endl;
                        }
                    }

                    string DFN_PT_csv_name = "PTPara";

                    CreateOrEmptyFile(DFNFileName + "/" + DFN_PT_csv_name +
                                      ".csv");

                    AddLineToFile("./" + DFNFileName + "/" + DFN_PT_csv_name +
                                      ".csv",
                                  "NumParticles," +
                                      std::to_string(NumberParticles) + ",\n");
                    AddLineToFile(
                        "./" + DFNFileName + "/" + DFN_PT_csv_name + ".csv",
                        "NumTimeSteps," + std::to_string(NumSteps) + ",\n");
                    AddLineToFile("./" + DFNFileName + "/" + DFN_PT_csv_name +
                                      ".csv",
                                  "MolecularDiffusion," +
                                      DoubleNumberToScientificNotationString(
                                          MolecularDiffusionCoefficient) +
                                      ",\n");
                    AddLineToFile(
                        "./" + DFNFileName + "/" + DFN_PT_csv_name + ".csv",
                        "DeltaT," +
                            DoubleNumberToScientificNotationString(DeltaT) +
                            ",\n");
                    AddLineToFile("./" + DFNFileName + "/" + DFN_PT_csv_name +
                                      ".csv",
                                  "InjectionMethod_initialcondition," +
                                      InjectionMethod_initialcondition + ",\n");
                    AddLineToFile("./" + DFNFileName + "/" + DFN_PT_csv_name +
                                      ".csv",
                                  "If_OutputAllPTInformationOrFPTCurve," +
                                      std::to_string(
                                          If_OutputAllPTInformationOrFPTCurve) +
                                      ",\n");
                    // SpacingOfControlPlanes
                    // IfOutputVarianceOfDisplacementsEachStep
                    // IfInjectAtCustomedPlane
                    // CustomedPlaneInjection
                    // IfUseFluxWeightedOrEqualProbableMixingIntersection
                    // TimeIntervalOutPTInformation
                    AddLineToFile(
                        "./" + DFNFileName + "/" + DFN_PT_csv_name + ".csv",
                        "SpacingOfControlPlanes," +
                            std::to_string(SpacingOfControlPlanes) + ",\n");
                    AddLineToFile(
                        "./" + DFNFileName + "/" + DFN_PT_csv_name + ".csv",
                        "IfOutputVarianceOfDisplacementsEachStep," +
                            std::to_string(
                                IfOutputVarianceOfDisplacementsEachStep) +
                            ",\n");
                    AddLineToFile(
                        "./" + DFNFileName + "/" + DFN_PT_csv_name + ".csv",
                        "IfInjectAtCustomedPlane," +
                            std::to_string(IfInjectAtCustomedPlane) + ",\n");
                    AddLineToFile(
                        "./" + DFNFileName + "/" + DFN_PT_csv_name + ".csv",
                        "CustomedPlaneInjection," +
                            std::to_string(CustomedPlaneInjection) + ",\n");
                    AddLineToFile(
                        "./" + DFNFileName + "/" + DFN_PT_csv_name + ".csv",
                        "IfUseFluxWeightedOrEqualProbableMixingIntersection," +
                            std::to_string(
                                IfUseFluxWeightedOrEqualProbableMixingIntersection) +
                            ",\n");
                    AddLineToFile(
                        "./" + DFNFileName + "/" + DFN_PT_csv_name + ".csv",
                        "TimeIntervalOutPTInformation," +
                            std::to_string(TimeIntervalOutPTInformation) +
                            ",\n");

                    AddLineToFile(
                        "./" + DFNFileName + "/" + DFN_PT_csv_name + ".csv",
                        "IfOutputAllParticlesAccumulativeDisplacements," +
                            std::to_string(
                                IfOutputAllParticlesAccumulativeDisplacements) +
                            ",\n");
                } // finish writting PTPara.csv

                cout << DFNFileName << ": PT running ...\n";
                //-- start run PT
                string DFN_PT_run_command = "cd ./" + DFNFileName + " && " +
                                            ExeuctablePath + "/DFN_PT " +
                                            "./PTPara";

                bool DFN_PT_success = RunCMD_with_RealTimeCheck(
                    DFN_PT_run_command, "./" + DFNFileName + "/" + LogFile);
                // cout << "DFN_PT_success: " << DFN_PT_success << endl;
                if (!DFN_PT_success)
                {
                    std::vector<uint> tmp_count(1);

                    if (!IfAFileExist("./" + DFNFileName + "/PT_ErrorCount.h5"))
                    {
                        cout << DFNFileName
                             << ": DFN_PT throw the first error, creating "
                                "`PT_ErrorCount.h5`\n";

                        h5g.NewFile("./" + DFNFileName + "/PT_ErrorCount.h5");

                        tmp_count[0] = 1;
                        h5g.AddDataset("./" + DFNFileName + "/PT_ErrorCount.h5",
                                       "N", "ErrorCount", tmp_count.data(),
                                       make_uint2(1, 0));
                    }
                    else
                    {
                        cout << DFNFileName
                             << ": DFN_PT throw an error again, reading "
                                "`PT_ErrorCount.h5`\n";
                        tmp_count = h5g.ReadDataset<uint>(
                            "./" + DFNFileName + "/PT_ErrorCount.h5", "N",
                            "ErrorCount");
                        tmp_count[0]++;
                        h5g.OverWrite("./" + DFNFileName + "/PT_ErrorCount.h5",
                                      "N", "ErrorCount", tmp_count.data(),
                                      make_uint2(1, 0));
                    }
                    cout
                        << DFNFileName
                        << ": DFN_PT failed. Found an exception. Error times = "
                        << tmp_count[0] << "\n";
                }
            }
        }
        catch (cuDFNsys::ExceptionsIgnore e)
        {
            cout << "cuDFNsys::Exceptions: " << e.what() << endl;
            result_system = system(
                ("rm -rf ./" + DFNFileName + "/*.h5 " + DFNFileName +
                 "/ParticlePositionResult/* ./" + DFNFileName + "/PTFinished")
                    .c_str());
            i--;
            continue;
        }
        catch (cuDFNsys::ExceptionsPause e)
        {
            cout << "cuDFNsys::Exceptions: " << e.what() << endl;
            result_system = system(
                ("rm -rf ./" + DFNFileName + "/*.h5 " + DFNFileName +
                 "/ParticlePositionResult/* ./" + DFNFileName + "/PTFinished")
                    .c_str());
            i--;
            continue;
        }
        catch (H5::Exception e)
        {
            cout << "H5::Exceptions: " << e.getDetailMsg() << endl;
            result_system = system(
                ("rm -rf ./" + DFNFileName + "/*.h5 " + DFNFileName +
                 "/ParticlePositionResult/* ./" + DFNFileName + "/PTFinished")
                    .c_str());
            i--;
            continue;
        }
        catch (H5::FileIException e)
        {
            cout << "H5::Exceptions: " << e.getDetailMsg() << endl;
            result_system = system(
                ("rm -rf ./" + DFNFileName + "/*.h5 " + DFNFileName +
                 "/ParticlePositionResult/* ./" + DFNFileName + "/PTFinished")
                    .c_str());
            i--;
            continue;
        }
        catch (...)
        {
            cout << "Unknown exceptions\n";
            result_system = system(
                ("rm -rf ./" + DFNFileName + "/*.h5 " + DFNFileName +
                 "/ParticlePositionResult/* ./" + DFNFileName + "/PTFinished")
                    .c_str());
            i--;
            continue;
        }
    }

    return 0;
}

bool RunCMD_with_RealTimeCheck(const string &cmd, const string &logFile,
                               const bool &IfGenLogFile)
{
    if (IfGenLogFile)
        CreateOrEmptyFile(logFile);

    FILE *pipe = popen((cmd + " 2>&1").c_str(), "r");
    if (!pipe)
    {
        //std::cerr << "popen() failed!";
        cout << "command failed: " << cmd << endl;
        return false;
    }

    // Read the output line by line
    char buffer[1024];
    bool statu = true;

    while (fgets(buffer, 1024, pipe) != nullptr)
    {
        std::string output(buffer);
        //std::cout << output; // Print the output in real-time
        AddLineToFile(logFile, output);
        // Check if the output contains any warning
        if (output.find("Warning") != std::string::npos)
        {
            //std::cerr << "Warning detected!\n";
            // Take appropriate action
            pclose(pipe);
            statu = false;
            return statu;
        }
    }

    // Close the pipe
    pclose(pipe);
    return statu;
}

void CreateOrEmptyFile(const std::string &filename)
{
    std::ofstream file(filename, std::ios::out | std::ios::trunc);
    if (file.is_open())
    {
        //std::cout << "File created or emptied successfully: " << filename
        //         << std::endl;
        file.close();
    }
    else
    {
        std::cerr << "Failed to create or empty file: " << filename
                  << std::endl;
    }
}
void AddLineToFile(const std::string &filename, const std::string &line)
{
    std::ofstream file(filename, std::ios::out | std::ios::app);
    if (file.is_open())
    {
        file << line;
        //std::cout << "Line added to file: " << filename << std::endl;
        file.close();
    }
    else
    {
        std::cerr << "Failed to open file: " << filename << std::endl;
    }
}
bool IfAFileExist(const string &FilePath)
{
    fs::path filePath = FilePath;

    // Check if the file exists
    if (fs::exists(filePath))
    {
        return true;
    }
    else
    {
        return false;
    }
}
string DoubleNumberToScientificNotationString(const double &number)
{
    std::stringstream ss;
    ss << std::scientific << number;
    std::string result = ss.str();
    return result;
}

uint GetH5DatasetSize(const string &nameH5, const string &nameDataset)
{
    H5File file(nameH5, H5F_ACC_RDONLY);

    // Open the dataset
    DataSet dataset = file.openDataSet(nameDataset);

    // Get the dataspace of the dataset
    DataSpace dataspace = dataset.getSpace();

    // Get the number of dimensions in the dataspace
    int ndims = dataspace.getSimpleExtentNdims();

    // Create a vector to store the size of each dimension
    hsize_t dims[ndims];

    // Get the size of each dimension
    dataspace.getSimpleExtentDims(dims);

    // Calculate the total size of the dataset
    uint totalSize = 1;
    for (int i = 0; i < ndims; ++i)
    {
        totalSize *= dims[i];
    }

    // Output the total size
    //std::cout << "Size of dataset: " << totalSize << std::endl;

    // Close the dataset and file
    dataset.close();
    file.close();

    return totalSize;
}
