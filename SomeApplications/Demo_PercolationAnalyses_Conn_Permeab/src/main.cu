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

void addLineToFile(const std::string &filename, const std::string &line);

int main(int argc, char *argv[])
{
    double i_start = cuDFNsys::CPUSecond();

    double domain_size = atof(argv[1]);
    int MCTimes = atoi(argv[2]);
    int Nproc = atoi(argv[3]);
    double WaitTimeCheckTasks = atof(argv[4]);
    if (WaitTimeCheckTasks < 1e-3)
    {
        cout << "WaitTimeCheckTasks cannot less than 1e-3. I will set it to be "
                "1e-3\n";
        WaitTimeCheckTasks = 1e-3;
    }

    std::vector<int2> TagLabel(SizeVec * MCTimes);

    for (int i = 0; i < SizeVec; ++i)
        for (int j = 0; j < MCTimes; ++j)
            TagLabel[i * MCTimes + j].x = i, TagLabel[i * MCTimes + j].y = j;

    std::random_device rd; // Obtain a random number from hardware
    std::default_random_engine eng(rd());                // Seed the generator
    std::shuffle(TagLabel.begin(), TagLabel.end(), eng); // Shuffle the vector

    string MainFile =
        "DFN_result_L_" + cuDFNsys::ToStringWithWidth(domain_size, 4);
    int resultas = system(string("mkdir " + MainFile).c_str());
    if (resultas != 0)
        cout << MainFile + " exists\n";

    //------------------
    string GetErrorFile = "./GetErrorFile.txt";
    string RecordProgressFile = "./RecordProgress.txt";
    system(string("touch " + GetErrorFile).c_str());
    system(string("echo '' > " + GetErrorFile).c_str());
    system(string("touch " + RecordProgressFile).c_str());
    system(string("echo '' > " + RecordProgressFile).c_str());

    double FinishedTasksCount = 0;
    std::vector<double> Bar_t(20);
    double start_value = 0.05, increment = 0.05;
    std::generate(Bar_t.begin(), Bar_t.end(),
                  [&start_value, increment]() mutable
                  {
                      double value = start_value;
                      start_value += increment;
                      return value;
                  });
    std::reverse(Bar_t.begin(), Bar_t.end());

#pragma omp parallel for schedule(dynamic) num_threads(Nproc)
    for (int i = 0; i < TagLabel.size(); ++i)
    {
        string DFN_File_path =
            "DFN_L_" + cuDFNsys::ToStringWithWidth(domain_size, 4) +
            "_RhoInd_" + cuDFNsys::ToStringWithWidth(TagLabel[i].x, 3) +
            "_MC_" + cuDFNsys::ToStringWithWidth(TagLabel[i].y, 5);

        system(string("mkdir -p " + MainFile + "/" + DFN_File_path).c_str());

        string Exe_run_command_string =
            " cd ./" + MainFile + "/" + DFN_File_path +
            " && bash ../../RunLevel2.sh " + DFN_File_path + " ../../Level2 " +
            (std::to_string(domain_size) + " ") +
            (std::to_string(TagLabel[i].x) + " ") +
            (std::to_string(TagLabel[i].y) + " ");

        std::system(Exe_run_command_string.c_str());

        int stop_idx = 0;

        while (stop_idx == 0)
        {

            std::this_thread::sleep_for(std::chrono::milliseconds(
                int(1000 * WaitTimeCheckTasks))); // 1000 milliseconds = 1 s

            //--------check log.txt
            if (fs::exists(string("./" + MainFile + "/" + DFN_File_path +
                                  "/Finished")))
            {
#pragma critical
                {
                    FinishedTasksCount++;
                    for (int k = 0; k < Bar_t.size(); ++k)
                        if (abs(FinishedTasksCount / (double(TagLabel.size())) -
                                Bar_t[k]) < 1e-6)
                        {

                            addLineToFile(
                                RecordProgressFile,
                                string(std::to_string(Bar_t[k] * 100) +
                                       "% have been finished\n"));
                            Bar_t.erase(Bar_t.begin() + k);
                            break;
                        }
                }
                stop_idx = 1;
                break;
            }

            //---------------------check if it is running by process name
            int result = system(string("pgrep -f " + DFN_File_path).c_str());
            if (result != 0)
            {
// record
#pragma critical
                {
                    addLineToFile(GetErrorFile,
                                  string(DFN_File_path +
                                         " is not running, checking by '" +
                                         string("pgrep -f " + DFN_File_path) +
                                         "'\n"));
                }
                // kill process by name
                system(string("pkill -f " + DFN_File_path).c_str());
                system(string("cd ./" + MainFile + "/" + DFN_File_path +
                              " && rm -rf ./Data.h5 ./Finished")
                           .c_str());
                stop_idx = 1;
                break;
            }

            //------------------check warning in log.txt
            std::string filename3 =
                "./" + MainFile + "/" + DFN_File_path + "/log.txt";
            {
                std::ifstream inputFile(filename3); // Open the file

                if (!inputFile)
                {
#pragma critical
                    {

                        addLineToFile(
                            GetErrorFile,
                            string("Unable to open file: " + filename3 + "," +
                                   Exe_run_command_string + "\n"));
                    }
                    // kill process by name
                    system(string("pkill -f " + DFN_File_path).c_str());
                    system(string("cd ./" + MainFile + "/" + DFN_File_path +
                                  " && rm -rf ./Data.h5 ./Finished")
                               .c_str());

                    stop_idx = 1;
                    inputFile.close(); // Close the file
                    break;
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

                                addLineToFile(
                                    GetErrorFile,
                                    string(
                                        "The log.txt contains the word "
                                        "'Warning', "
                                        "I am going to delete the hdf5 and the "
                                        "Finished" +
                                        Exe_run_command_string + "\n"));
                            }
                            // kill process by name
                            system(string("pkill -f " + DFN_File_path).c_str());

                            system(string("cd ./" + MainFile + "/" +
                                          DFN_File_path +
                                          " && rm -rf ./Data.h5 ./Finished")
                                       .c_str());

                            stop_idx = 1;
                            break; // Exit loop once the word is found
                        }
                    }
                }

                inputFile.close(); // Close the file
            }
        };
    }

    double ElapseTime = cuDFNsys::CPUSecond() - i_start;

    cuDFNsys::HDF5API hjk5;

    // Get the current time point
    std::chrono::system_clock::time_point now =
        std::chrono::system_clock::now();

    // Convert the time point to a time_t object
    std::time_t now_c = std::chrono::system_clock::to_time_t(now);

    // Convert the time_t object to a string representation
    std::string timeString = std::ctime(&now_c);
    for (char &c : timeString)
    {
        if (c == ' ')
        {
            c = '_';
        }
    }

    std::string filenameERD =
        "TotalRunTime_L_" + cuDFNsys::ToStringWithWidth(domain_size, 4) + ".h5";
    if (!fs::exists(filenameERD))
        hjk5.NewFile(filenameERD);

    hjk5.AddDataset(filenameERD, "N", "TotalRunTime_" + timeString, &ElapseTime,
                    make_uint2(1, 0));

    return 0;
};

void addLineToFile(const std::string &filename, const std::string &line)
{
    // Open the file in append mode
    std::ofstream file(filename, std::ios::app);

    if (file.is_open())
    {
        // Write a line to the file
        file << line;

        // Close the file
        file.close();

        //std::cout << "Line added to file." << std::endl;
    }
    else
    {
        //std::cerr << "Unable to open file." << std::endl;
        // Consider throwing an exception or handling the error in some other way
    }
}