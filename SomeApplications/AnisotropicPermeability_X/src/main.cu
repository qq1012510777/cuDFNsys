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
    int PercoDir = 0;
    double Kappa = atof(argv[6]);

    int FracNumInit = atoi(argv[7]);
    int FracNumIncre = atoi(argv[8]);
    int NumFracIncre = atoi(argv[9]);

    int FracSizeDistriMode = atoi(argv[10]);
    double4 FracSizeDistriPara = make_double4(atof(argv[11]), atof(argv[12]),
                                              atof(argv[13]), atof(argv[14]));

    double Beta = atof(argv[15]);
    double Gamma = atof(argv[16]);
    // DFN mesh parameters
    double MeshMinimumGridSize = atof(argv[17]);
    double MeshMaximumGridSize = atof(argv[18]);

    // exit(0);
    cuDFNsys::HDF5API h5g;

    for (int i = 0; i <= NumFracIncre; ++i)
    {

        string DFNFileName = "DFN_X_" + cuDFNsys::ToStringWithWidth(i + 1, 3);

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

                    if (!DFN_gen_success)
                    {
                        cout << DFNFileName << ": DFN_Gen failed\n";
                        result_system =
                            system(("rm -rf ./" + DFNFileName + "/Class_DFN.h5")
                                       .c_str());
                        continue;
                    }

                    system(("cp ./" + DFNFileName + "/Class_DFN.h5  ./" +
                            DFNFileName + "/Class_DFN_original.h5")
                               .c_str());
                    break;
                }
            }

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

            if (GetH5DatasetSize("./" + DFNFileName + "/Class_DFN.h5",
                                 "PercolationClusters") > 0)
            {
                //-----------------DFN_Mesh------------
                if (IfAFileExist("./" + DFNFileName + "/Class_DFN.h5") &&
                    !IfAFileExist("./" + DFNFileName + "/Class_MESH.h5"))
                {

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

                //------------------
            }
            CreateOrEmptyFile("./" + DFNFileName + "/X_DarcyFlow_Finished");
        }
        catch (cuDFNsys::ExceptionsIgnore e)
        {
            cout << "cuDFNsys::Exceptions: " << e.what() << endl;
            result_system = system(
                ("rm -rf ./" + DFNFileName + "/*.h5 " + DFNFileName +
                 "/ParticlePositionResult/* ./" + DFNFileName + "/X_DarcyFlow_Finished")
                    .c_str());
            i--;
            continue;
        }
        catch (cuDFNsys::ExceptionsPause e)
        {
            cout << "cuDFNsys::Exceptions: " << e.what() << endl;
            result_system = system(
                ("rm -rf ./" + DFNFileName + "/*.h5 " + DFNFileName +
                 "/ParticlePositionResult/* ./" + DFNFileName + "/X_DarcyFlow_Finished")
                    .c_str());
            i--;
            continue;
        }
        catch (H5::Exception e)
        {
            cout << "H5::Exceptions: " << e.getDetailMsg() << endl;
            result_system = system(
                ("rm -rf ./" + DFNFileName + "/*.h5 " + DFNFileName +
                 "/ParticlePositionResult/* ./" + DFNFileName + "/X_DarcyFlow_Finished")
                    .c_str());
            i--;
            continue;
        }
        catch (H5::FileIException e)
        {
            cout << "H5::Exceptions: " << e.getDetailMsg() << endl;
            result_system = system(
                ("rm -rf ./" + DFNFileName + "/*.h5 " + DFNFileName +
                 "/ParticlePositionResult/* ./" + DFNFileName + "/X_DarcyFlow_Finished")
                    .c_str());
            i--;
            continue;
        }
        catch (...)
        {
            cout << "Unknown exceptions\n";
            result_system = system(
                ("rm -rf ./" + DFNFileName + "/*.h5 " + DFNFileName +
                 "/ParticlePositionResult/* ./" + DFNFileName + "/X_DarcyFlow_Finished")
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
        // std::cerr << "popen() failed!";
        cout << "command failed: " << cmd << endl;
        return false;
    }

    // Read the output line by line
    char buffer[1024];
    bool statu = true;

    while (fgets(buffer, 1024, pipe) != nullptr)
    {
        std::string output(buffer);
        // std::cout << output; // Print the output in real-time
        AddLineToFile(logFile, output);
        // Check if the output contains any warning
        if (output.find("Warning") != std::string::npos)
        {
            // std::cerr << "Warning detected!\n";
            //  Take appropriate action
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
        // std::cout << "File created or emptied successfully: " << filename
        //          << std::endl;
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
        // std::cout << "Line added to file: " << filename << std::endl;
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
    // std::cout << "Size of dataset: " << totalSize << std::endl;

    // Close the dataset and file
    dataset.close();
    file.close();

    return totalSize;
}
