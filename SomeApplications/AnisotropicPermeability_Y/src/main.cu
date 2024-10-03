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

    int NumFracIncre = atoi(argv[2]);

    // exit(0);
    cuDFNsys::HDF5API h5g;

    for (int i = 0; i <= NumFracIncre; ++i)
    {

        string DFNFileName = "DFN_Y_" + cuDFNsys::ToStringWithWidth(i + 1, 3);

        string DFNFileName_X = "DFN_X_" + cuDFNsys::ToStringWithWidth(i + 1, 3);

        int result_system = system(("mkdir ./" + DFNFileName + " -p").c_str());

        string LogFile =
            "log_DFN_" + cuDFNsys::ToStringWithWidth(i + 1, 3) + ".txt";

        CreateOrEmptyFile("./" + DFNFileName + "/" + LogFile);
        try
        {
            if (IfAFileExist("./" + DFNFileName + "/Y_DarcyFlow_Finished"))
                continue;

            //-----------------DFN_Gen-------------
            if (IfAFileExist("./" + DFNFileName_X + "/Class_DFN_original.h5"))
            {

                result_system = system(("cp -f ./" + DFNFileName_X + "/Class_DFN_original.h5   ./" + DFNFileName + "/Class_DFN_original.h5")
                                           .c_str());

                string DFN_gen_run_command = "cd ./" + DFNFileName + " && " +
                                             ExeuctablePath + "/LoadDFNAndChangePercolationDirection " +
                                             "./Class_DFN_original  1";

                bool DFN_gen_success = RunCMD_with_RealTimeCheck(
                    DFN_gen_run_command, "./" + DFNFileName + "/" + LogFile,
                    true);
                // cout << DFN_gen_run_command << endl;

                if (!DFN_gen_success)
                {
                    cout << DFNFileName << ": DFN_load failed\n";
                    result_system =
                        system(("rm -rf ./" + DFNFileName + "/Class_DFN.h5")
                                   .c_str());
                    continue;
                }
            }
            else
            {
                cout << DFNFileName << ": DFN_X does not have a DFN\n";
                continue;
            }

            std::vector<uint> PercoCluID_X = h5g.ReadDataset<uint>("./" + DFNFileName_X + "/Class_DFN.h5", "N", "PercolationClusters");
            std::vector<uint> PercoFracID_X;
            for (int j = 0; j < PercoCluID_X.size(); ++j)
            {
                std::vector<uint> tmp = h5g.ReadDataset<uint>("./" + DFNFileName_X + "/Class_DFN.h5", "N", "Cluster_" + std::to_string(PercoCluID_X[j] + 1));
                PercoFracID_X.insert(PercoFracID_X.end(), tmp.begin(), tmp.end());
            }

            if (GetH5DatasetSize("./" + DFNFileName + "/Class_DFN.h5",
                                 "PercolationClusters") > 0)
            {
                std::vector<uint> PercoCluID_Y = h5g.ReadDataset<uint>("./" + DFNFileName + "/Class_DFN.h5", "N", "PercolationClusters");
                std::vector<uint> PercoFracID_Y;
                for (int j = 0; j < PercoCluID_Y.size(); ++j)
                {
                    std::vector<uint> tmp = h5g.ReadDataset<uint>("./" + DFNFileName + "/Class_DFN.h5", "N", "Cluster_" + std::to_string(PercoCluID_Y[j] + 1));
                    PercoFracID_Y.insert(PercoFracID_Y.end(), tmp.begin(), tmp.end());
                }

                bool sameDFN = false;
                if (PercoFracID_X.size() == PercoFracID_Y.size() && std::equal(PercoFracID_X.begin(), PercoFracID_X.end(), PercoFracID_Y.begin()))
                {
                    sameDFN = true;
                    cout << "The DFN is the same with the one in the X direction\n";
                }
                else
                    cout << "The DFN is NOT the same with the one in the X direction\n";

                //-----------------DFN_Mesh------------
                if ((IfAFileExist("./" + DFNFileName + "/Class_DFN.h5") &&
                     !IfAFileExist("./" + DFNFileName + "/Class_MESH.h5") &&
                     !IfAFileExist("./" + DFNFileName_X + "/Class_MESH.h5")) ||
                    !sameDFN)
                {
                GenMesh:;
                    result_system = system(("cp -f ./" + DFNFileName_X + "/MeshPara.csv   ./" + DFNFileName + "/MeshPara.csv")
                                               .c_str());

                    string DFN_mesh_run_command = "cd ./" + DFNFileName + " && " +
                                                  ExeuctablePath + "/DFN_Mesh " +
                                                  "./MeshPara";
                    bool DFN_mesh_success = RunCMD_with_RealTimeCheck(
                        DFN_mesh_run_command, "./" + DFNFileName + "/" + LogFile);

                    if (!DFN_mesh_success)
                    {

                        cout << DFNFileName
                             << ": DFN_Mesh failed \n";
                        result_system =
                            system(("rm -rf ./" + DFNFileName + "/Class_DFN.h5 ./" +
                                    DFNFileName + "/Class_MESH.h5")
                                       .c_str());
                        continue;
                    }
                }
                else if (IfAFileExist("./" + DFNFileName + "/Class_DFN.h5") &&
                         !IfAFileExist("./" + DFNFileName + "/Class_MESH.h5") &&
                         IfAFileExist("./" + DFNFileName_X + "/Class_MESH.h5"))
                {
                    if (!sameDFN)
                        goto GenMesh;
                    result_system = system(("cp -f ./" + DFNFileName_X + "/Class_MESH.h5   ./" + DFNFileName + "/Class_MESH.h5")
                                               .c_str());

                    string DFN_mesh_load_command = "cd ./" + DFNFileName + " && " +
                                                   ExeuctablePath + "/LoadMeshRenumberingEdge  1";
                    bool DFN_mesh_success = RunCMD_with_RealTimeCheck(
                        DFN_mesh_load_command, "./" + DFNFileName + "/" + LogFile);
                    // cout << DFN_mesh_load_command << endl;
                    if (!DFN_mesh_success)
                    {

                        cout << DFNFileName
                             << ": DFN_Mesh_load failed\n";
                        result_system =
                            system(("rm -rf ./" + DFNFileName + "/Class_DFN.h5 ./" +
                                    DFNFileName + "/Class_MESH.h5")
                                       .c_str());
                        continue;
                    }
                }

                //-----------------DFN_Flow---------------
                if (IfAFileExist("./" + DFNFileName + "/Class_DFN.h5") &&
                    IfAFileExist("./" + DFNFileName + "/Class_MESH.h5") &&
                    !IfAFileExist("./" + DFNFileName + "/Class_FLOW.h5"))
                {
                    result_system = system(("cp -f ./" + DFNFileName_X + "/FlowPara.csv   ./" + DFNFileName + "/FlowPara.csv")
                                               .c_str());

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
            else
            {
                CreateOrEmptyFile("./" + DFNFileName + "/no_cluster_in_Y");
            }

            string DFN_SToreData_command = "cd ./" + DFNFileName + " && " +
                                           "python3 " + ExeuctablePath + "/StoreAndDelteHDF5.py " +
                                           "1";
            bool DFN_SToreData_command_success = RunCMD_with_RealTimeCheck(
                DFN_SToreData_command, "./" + DFNFileName + "/" + LogFile);
            if (!DFN_SToreData_command_success)
            {
                if (!DFN_SToreData_command_success)
                    cout << DFNFileName
                         << ": DFN store data failed, regenerate DFN\n";
                result_system =
                    system(("rm -rf ./" + DFNFileName + "/Class_DFN.h5 ./" +
                            DFNFileName + "/Class_MESH.h5 ./" +
                            DFNFileName + "/Class_FLOW.h5")
                               .c_str());
                i--;
                continue;
            }
            else
                result_system = system(
                    ("cd ./" + DFNFileName + " && rm -rf Class_FLOW.h5 Class_MESH.h5 DFN_FLOW_VISUAL.h5 DFN_MESH_VISUAL.h5 DFN_VISUAL.h5").c_str());

            CreateOrEmptyFile("./" + DFNFileName + "/Y_DarcyFlow_Finished");
        }
        catch (cuDFNsys::ExceptionsIgnore e)
        {
            cout << "cuDFNsys::Exceptions: " << e.what() << endl;
            result_system = system(
                ("rm -rf ./" + DFNFileName + "/*.h5 " + DFNFileName +
                 "/ParticlePositionResult/* ./" + DFNFileName + "/Y_DarcyFlow_Finished")
                    .c_str());
            i--;
            continue;
        }
        catch (cuDFNsys::ExceptionsPause e)
        {
            cout << "cuDFNsys::Exceptions: " << e.what() << endl;
            result_system = system(
                ("rm -rf ./" + DFNFileName + "/*.h5 " + DFNFileName +
                 "/ParticlePositionResult/* ./" + DFNFileName + "/Y_DarcyFlow_Finished")
                    .c_str());
            i--;
            continue;
        }
        catch (H5::Exception e)
        {
            cout << "H5::Exceptions: " << e.getDetailMsg() << endl;
            result_system = system(
                ("rm -rf ./" + DFNFileName + "/*.h5 " + DFNFileName +
                 "/ParticlePositionResult/* ./" + DFNFileName + "/Y_DarcyFlow_Finished")
                    .c_str());
            i--;
            continue;
        }
        catch (H5::FileIException e)
        {
            cout << "H5::Exceptions: " << e.getDetailMsg() << endl;
            result_system = system(
                ("rm -rf ./" + DFNFileName + "/*.h5 " + DFNFileName +
                 "/ParticlePositionResult/* ./" + DFNFileName + "/Y_DarcyFlow_Finished")
                    .c_str());
            i--;
            continue;
        }
        catch (...)
        {
            cout << "Unknown exceptions\n";
            result_system = system(
                ("rm -rf ./" + DFNFileName + "/*.h5 " + DFNFileName +
                 "/ParticlePositionResult/* ./" + DFNFileName + "/Y_DarcyFlow_Finished")
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
