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
#include <omp.h>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <unistd.h>
#include <vector>
using namespace std;

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
    if (argc < 4)
    {
        cout << "Please input the following arguments:\n";
        cout << "Argument 1: the inital set number\n";
        cout << "Argument 2: the end set number \n";
        cout << "Argument 3: the number of threads for parallelization\n";
    }

    int a = atoi(argv[1]);
    int b = atoi(argv[2]);
    int Nproc = atoi(argv[3]);

    cout << "Argument 1: " << a << "\n";
    cout << "Argument 2: " << b << "\n";
    cout << "Argument 3: " << Nproc << "\n";

    std::vector<int> SleepSeconds(Nproc, 0);
    std::iota(SleepSeconds.begin(), SleepSeconds.end(), 0);
    std::for_each(SleepSeconds.begin(), SleepSeconds.end(),
                  [](int &num) { num *= 10; });

#pragma omp parallel for schedule(dynamic) num_threads(Nproc)
    for (int i = a; i <= b; ++i)
    {
        int timeSleep = SleepSeconds[(i - a) % 3];

        std::this_thread::sleep_for(std::chrono::seconds(timeSleep));

        string command = "bash run.sh " + std::to_string(i) + " > log_set" +
                         addZeroPadding(i, 5) + ".txt 2>&1";

        cout << "Running " << command << " after sleep " << timeSleep
             << " seconds" << endl;
        GetTime();

        int result = system(command.c_str());

        // Check if the command executed successfully
        if (result == 0)
        {
            // Command executed successfully
            cout << command << " executed successfully\n";
        }
        else
        {
            // Command failed
            cout << command << " failed\n";
        }
    }
    return 0;
}