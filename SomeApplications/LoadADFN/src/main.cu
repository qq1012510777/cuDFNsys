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


int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        std::cout << "Usage: " << argv[0] << " <DFN H5 filename without the suffix>\n";
        exit(0);
    }
    cuDFNsys::DFN<double> my_dfn;
    my_dfn.LoadClassFromH5(std::string(argv[1]));

    std::string DFNVisualFIle = "DFN_Visual";
    my_dfn.Visualization(DFNVisualFIle, DFNVisualFIle, DFNVisualFIle, false, false, true, true);
    return 0;
}
