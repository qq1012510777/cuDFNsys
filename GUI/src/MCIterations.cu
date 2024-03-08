#include "cuDFNsys.cuh"
#include <algorithm>
#include <chrono>
#include <cstdlib> // For std::system
#include <ctime>
#include <fstream>
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

bool RunCMD_with_RealTimeCheck(const string &cmd, const string &logFile, const bool &IfGenLogFile = false);
void CreateOrEmptyFile(const std::string &filename);
void AddLineToFile(const std::string &filename, const std::string &line);
bool IfAFileExist(const string &FilePath);
void CriticalOutput(const string &Sstring)
{
#pragma critical
    {
        cout << Sstring;
    }
}

int main(int argc, char *argv[])
{
    string cuDFNsys_GUI_root = (argc >= 2 ? argv[1] : "~/cuDFNsys/GUI");
    int NumIterations = (argc >= 3 ? atoi(argv[2]) : 1);
    int NumProcessors = (argc >= 4 ? atoi(argv[3]) : 1);

    // ---------------generate a vector of random numbers -------------
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(1,
                                        5); // Generate integers between 1 and 5
    // Create a vector to store random numbers
    std::vector<int> randomNumbers(NumIterations);
    // Generate and store random numbers in the vector
    for (int i = 0; i < NumIterations; ++i)
        randomNumbers[i] = dis(gen);

    //-----------------------------------------------------------------
    const int PercentageIncrement = 5;
    int lastOutputPercentage = 0;
#pragma omp parallel for schedule(dynamic) num_threads(NumProcessors)
    for (int i = 0; i < NumIterations; ++i)
    {
        string ProjectName = "./MC_" + cuDFNsys::ToStringWithWidth(i + 1, 6);
        system(("mkdir " + ProjectName + " -p").c_str());
        string LogFile =
            "log_MC_" + cuDFNsys::ToStringWithWidth(i + 1, 6) + ".txt";

        if (IfAFileExist(ProjectName + "/Class_DFN.h5") &&
            IfAFileExist(ProjectName + "/Class_MESH.h5") &&
            IfAFileExist(ProjectName + "/Class_FLOW.h5"))
            continue;

        //--------------------DFN gen
        string DFN_gen_run_command = "cd ./" + ProjectName + " && " +
                                     cuDFNsys_GUI_root + "/DFN_Gen " +
                                     "../StochasticDFN";

        std::this_thread::sleep_for(std::chrono::seconds(
            randomNumbers[i])); // sleep a little bit for diverse random seeds

        bool DFN_gen_success =
            RunCMD_with_RealTimeCheck(DFN_gen_run_command, LogFile, true);

        if (!DFN_gen_success)
        {
            CriticalOutput("MC_" + cuDFNsys::ToStringWithWidth(i + 1, 6) +
                           " fails in DFN_Gen\n");
            continue;
        }

        //--------------------DFN mesh
        string DFN_mesh_run_command = "cd " + ProjectName + " && " +
                                      cuDFNsys_GUI_root + "/DFN_Mesh " +
                                      "../MeshPara";
        bool DFN_mesh_success =
            RunCMD_with_RealTimeCheck(DFN_mesh_run_command, LogFile);

        if (!DFN_mesh_success)
        {
            CriticalOutput("MC_" + cuDFNsys::ToStringWithWidth(i + 1, 6) +
                           " fails in DFN_Mesh\n");
            continue;
        }

        //--------------------DFN Flow
        string DFN_flow_run_command = "cd " + ProjectName + " && " +
                                      cuDFNsys_GUI_root + "/DFN_Flow " +
                                      "../FlowPara";
        bool DFN_flow_success =
            RunCMD_with_RealTimeCheck(DFN_flow_run_command, LogFile);

        if (!DFN_flow_success)
        {
            CriticalOutput("MC_" + cuDFNsys::ToStringWithWidth(i + 1, 6) +
                           " fails in DFN_Flow\n");
            continue;
        }
#pragma critical
        {
            // Calculate progress percentage
            double progress =
                (static_cast<double>(i + 1) / NumIterations) * 100;

            // Check if the progress matches the next equidistant percentage increment
            if (static_cast<int>(progress) >=
                lastOutputPercentage + PercentageIncrement)
            {
                lastOutputPercentage += PercentageIncrement;
                std::cout << "Progress: " << lastOutputPercentage
                          << "% completed" << std::endl;
            }
        }
    }

    cout << "\nAll MC finished\n";

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
        return false;
    }

    // Read the output line by line
    char buffer[1024];
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
            return false;
        }
    }

    // Close the pipe
    pclose(pipe);
    return true;
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