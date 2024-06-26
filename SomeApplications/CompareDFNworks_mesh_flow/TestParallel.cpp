#include <algorithm>
#include <chrono>
#include <cstdlib> // For std::system
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits.h>
#include <memory>
#include <numeric>
#include <omp.h>
#include <random>
#include <sstream>
#include <string>
#include <sys/time.h>
#include <thread>
#include <unistd.h>
#include <vector>
using namespace std;

double CPUSecond()
{
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return ((double)tp.tv_sec + (double)tp.tv_usec * 1e-6);
}; // CPUSecond

std::string addZeroPadding(int number, int padding)
{
    // Convert number to string
    std::ostringstream oss;
    oss << number;
    std::string numStr = oss.str();

    // Calculate padding length
    int numDigits = numStr.length();
    int paddingLength = padding - numDigits;

    // Add zero padding if needed
    if (paddingLength > 0)
    {
        numStr = std::string(paddingLength, '0') + numStr;
    }

    return numStr;
}

void GetTime()
{
    // Get the current time_point
    std::chrono::system_clock::time_point now =
        std::chrono::system_clock::now();

    // Convert time_point to time_t
    std::time_t now_c = std::chrono::system_clock::to_time_t(now);

    // Convert time_t to local time
    std::tm *local_time = std::localtime(&now_c);

    // Print the current date and time
    std::cout << "Current date and time: ";
    std::cout << local_time->tm_year + 1900 << '-' << local_time->tm_mon + 1
              << '-' << local_time->tm_mday << ' ' << local_time->tm_hour << ':'
              << local_time->tm_min << ':' << local_time->tm_sec << std::endl;
}

int main(int argc, char *argv[])
{
    std::filesystem::path current_path_oo = std::filesystem::current_path();

    std::string current_path = current_path_oo.string();

    for (int i = 1; i <= 10; ++i)
    {
        string threadFile =
            current_path + "/ThreadParallel_" + std::to_string(i);

        int resultsq = system(("mkdir -p " + threadFile).c_str());

        double time_start = CPUSecond();

        cout << "using CPU " << i << " threads in parallel\n";
#pragma omp parallel for schedule(dynamic) num_threads(i)
        for (int j = 0; j < 30; ++j)
        {
            string command =
                "cd " + threadFile + " && " + current_path + "/main " +
                std::to_string(j + 1) + " " + current_path +
                "/example_FracPara_PercolativeFractures  > ./log_" +
                std::to_string(j + 1) + ".txt 2>&1";

            int result = system(command.c_str());

            // Check if the command executed successfully
            if (result == 0)
            {
                // Command executed successfully
                //cout << command << " executed successfully\n";
            }
            else
            {
                // Command failed
                cout << command << " failed\n";
                exit(0);
            }
        }
        double time_elapse = CPUSecond() - time_start;

        string timeFile = threadFile + "/TotalTime.txt";

        // Create and open a text file
        std::ofstream outfile(timeFile.c_str());

        // Check if the file is open
        if (outfile.is_open())
        {
            // Write the number to the file
            outfile << time_elapse;

            // Close the file
            outfile.close();
        }
        else
        {
            std::cerr << "Error opening file." << timeFile << std::endl;
        }

        cout << "Running time: " << time_elapse << " seconds\n";
    }
    return 0;
}