#include "cuDFNsys.cuh"
#include <algorithm>
#include <chrono>
#include <cstdlib> // For std::system
#include <ctime>
#include <fstream>
#include <iostream>
#include <limits.h>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <unistd.h>
#include <vector>

namespace fs = std::filesystem;

using namespace std;

const int SizeVec = 30;

int main(int argc, char *argv[])
{
    double i_start = cuDFNsys::CPUSecond();

    double domain_size = atof(argv[1]);
    int MCTimes = atoi(argv[2]);

    std::vector<int2> TagLabel(SizeVec * MCTimes);

    for (int i = 0; i < SizeVec; ++i)
        for (int j = 0; j < MCTimes; ++j)
            TagLabel[i * MCTimes + j].x = i, TagLabel[i * MCTimes + j].y = j;

    std::random_device rd; // Obtain a random number from hardware
    std::default_random_engine eng(rd());                // Seed the generator
    std::shuffle(TagLabel.begin(), TagLabel.end(), eng); // Shuffle the vector

    string MainFile =
        "DFN_result_L_" + cuDFNsys::ToStringWithWidth(domain_size, 4);

    string command1 = "mkdir " + MainFile;
    int resultas = system(command1.c_str());

    if (resultas != 0)
    {
        cout << MainFile + " exists\n";
    }

    //------------------
    string GetErrorFile = "./GetErrorFile.txt";
    string command_78 = "touch " + GetErrorFile;
    system(command_78.c_str());
    command_78 = "echo '' > " + GetErrorFile;
    system(command_78.c_str());

    std::ofstream outputFile(GetErrorFile, std::ios::app);
    if (!outputFile.is_open())
    {
        string command_Qd =
            "echo CANNOT_OPEN_GetErrorFile_TEXT > " + GetErrorFile;
        system(command_Qd.c_str());
        outputFile.close();
        return 0;
    }

#pragma omp parallel for schedule(dynamic) num_threads(8)
    for (int i = 0; i < TagLabel.size(); ++i)
    {
        string path2 = "DFN_L_" + cuDFNsys::ToStringWithWidth(domain_size, 4) +
                       "_RhoInd_" +
                       cuDFNsys::ToStringWithWidth(TagLabel[i].x, 3) + "_MC_" +
                       cuDFNsys::ToStringWithWidth(TagLabel[i].y, 5);
        string command1 = "mkdir -p " + MainFile + "/" + path2;
        system(command1.c_str());

        string Input = " cd ./" + MainFile + "/" + path2 +
                       " && bash ../../RunLevel2.sh " + path2 +
                       " ../../Level2 ";
        Input += (std::to_string(domain_size) + " ");
        Input += (std::to_string(TagLabel[i].x) + " ");
        Input += (std::to_string(TagLabel[i].y) + " ");
        //Input += (" ");

        // Get the current time point
        std::chrono::system_clock::time_point now =
            std::chrono::system_clock::now();

        // Convert the time point to a time_t object
        std::time_t now_c = std::chrono::system_clock::to_time_t(now);

        // Convert the time_t object to a string representation
        std::string timeString = std::ctime(&now_c);

#pragma critical
        {
            cout
                << "Running [" << Input << "] at time: " << timeString
                << endl; // this just a output, so this would not take a lot of time
        }
        std::system(Input.c_str());

        int stop_idx = 0;

        while (stop_idx == 0)
        {
            std::this_thread::sleep_for(std::chrono::seconds(0.8));

            //--------check log.txt
            std::string filename2 = "./" + MainFile + "/" + path2 + "/Finished";
            if (fs::exists(filename2))
            {
#pragma critical
                {
                    std::cout << "[" << Input << "] finished!" << std::endl;
                }
                stop_idx = 1;
            }

            //--------check data.h5
            std::string filename16 = "./" + MainFile + "/" + path2 + "/Data.h5";
            if (!fs::exists(filename16))
            {
#pragma critical
                {
                    //std::cout << Input << " was not running!\n" << std::endl;
                    string command_qw = Input + " was not running!\n";
                    outputFile << command_qw;
                }
                stop_idx = 1;
            }

            //------------------check warning in log.txt
            std::string filename3 = "./" + MainFile + "/" + path2 + "/log.txt";
            {
                std::ifstream inputFile(filename3); // Open the file

                if (!inputFile)
                {
#pragma critical
                    {
                        outputFile << "Unable to open file: " << filename3
                                   << "," << Input << endl;
                    }
                    stop_idx = 1;
                }
                else
                {

                    std::string line;
                    while (std::getline(inputFile, line))
                    { // Read the file line by line
                        // Search for the word "warning" in each line
                        if (line.find("Warning") != std::string::npos)
                        {
#pragma critical
                            {
                                outputFile
                                    << "The log.txt contains the word 'Warning'"
                                    << ", I am going to delete the hdf5 and "
                                       "the Finished"
                                    << Input << endl;
                            }
                            string command_3 = "pkill -f " + path2;
                            system(command_3.c_str());

                            string command_4 =
                                "cd ./" + MainFile + "/" + path2 +
                                " && gio trash -f Data.h5 Finished";
                            system(command_4.c_str());
                            stop_idx = 1;
                            break; // Exit loop once the word is found
                        }
                    }
                }

                inputFile.close(); // Close the file
            }
        };
        outputFile.close();
    }

    double ElapseTime = cuDFNsys::CPUSecond() - i_start;

    cuDFNsys::HDF5API hjk5;
    hjk5.NewFile("TotalRunTime_L_" +
                 cuDFNsys::ToStringWithWidth(domain_size, 4) + ".h5");
    hjk5.AddDataset("TotalRunTime_L_" +
                        cuDFNsys::ToStringWithWidth(domain_size, 4) + ".h5",
                    "N", "TotalRunTime", &ElapseTime, make_uint2(1, 0));

    return 0;
};