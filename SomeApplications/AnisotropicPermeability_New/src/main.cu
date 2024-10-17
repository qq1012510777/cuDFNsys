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

string ExeuctablePath;
// DFN generation parameters
double DomainSizeX;
double3 DomainDimensionRatio;

double Kappa;
int FracNumInit;
int FracNumIncre;
int NumFracIncre;
int FracSizeDistriMode;
double4 FracSizeDistriPara;
double Beta;
double Gamma;
// DFN mesh parameters
double MeshMinimumGridSize;
double MeshMaximumGridSize;

void GenCsvFiles();

int main(int argc, char *argv[])
{

    string cwd = fsm::current_path().string();

    time_t t;
    time(&t);

    ExeuctablePath = argv[1];
    // DFN generation parameters
    DomainSizeX = atof(argv[2]);
    DomainDimensionRatio =
        make_double3(atof(argv[3]), atof(argv[4]), atof(argv[5]));
    PercoDir = 0;
    Kappa = atof(argv[6]);

    FracNumInit = atoi(argv[7]);
    FracNumIncre = atoi(argv[8]);
    NumFracIncre = atoi(argv[9]);

    FracSizeDistriMode = atoi(argv[10]);
    FracSizeDistriPara = make_double4(atof(argv[11]), atof(argv[12]),
                                      atof(argv[13]), atof(argv[14]));

    Beta = atof(argv[15]);
    Gamma = atof(argv[16]);
    // DFN mesh parameters
    MeshMinimumGridSize = atof(argv[17]);
    MeshMaximumGridSize = atof(argv[18]);

    // exit(0);
    cuDFNsys::HDF5API h5g;

    for (int i = 0; i <= NumFracIncre; ++i)
    {
        string DFNFileName = cwd + "/DFN_" + cuDFNsys::ToStringWithWidth(i + 1, 3);
        int result_system = system(("mkdir " + DFNFileName + " -p").c_str());

        for (int j = 0; j < 3; ++j)
        {

        }
    }

    return 0;
}

void GenCsvFiles(const string &DFNFileName, const int &stepNo_i)
{
    string DFN_Gen_csv_name = "StochasticDFN";

    CreateOrEmptyFile(DFNFileName + "/" + DFN_Gen_csv_name +
                      ".csv");

    AddLineToFile("" + DFNFileName + "/" + DFN_Gen_csv_name +
                      ".csv",
                  "IfStochastic,1,\n");
    AddLineToFile(
        "" + DFNFileName + "/" + DFN_Gen_csv_name + ".csv",
        "DomainSizeX," + std::to_string(DomainSizeX) + ",\n");
    AddLineToFile(
        "" + DFNFileName + "/" + DFN_Gen_csv_name + ".csv",
        "DomainDimensionRatio," +
            std::to_string(DomainDimensionRatio.x) + "," +
            std::to_string(DomainDimensionRatio.y) + "," +
            std::to_string(DomainDimensionRatio.z) + ",\n");
    AddLineToFile("" + DFNFileName + "/" + DFN_Gen_csv_name +
                      ".csv",
                  "Percolation_direction, 0,\n");
    AddLineToFile("" + DFNFileName + "/" + DFN_Gen_csv_name +
                      ".csv",
                  "NumFractureGroups,1,\n");
    AddLineToFile(
        "" + DFNFileName + "/" + DFN_Gen_csv_name + ".csv",
        "NumFractureEachGroup," +
            std::to_string(FracNumInit + stepNo_i * FracNumIncre) + ",\n");
    AddLineToFile("" + DFNFileName + "/" + DFN_Gen_csv_name +
                      ".csv",
                  "KappaValues," + std::to_string(Kappa) + ",\n");
    AddLineToFile("" + DFNFileName + "/" + DFN_Gen_csv_name +
                      ".csv",
                  "MeanOrientationOfFisherDistribution,0,0,1,\n");
    AddLineToFile("" + DFNFileName + "/" + DFN_Gen_csv_name +
                      ".csv",
                  "ModeOfSizeDistribution," +
                      std::to_string(FracSizeDistriMode) + ",\n");
    AddLineToFile("" + DFNFileName + "/" + DFN_Gen_csv_name +
                      ".csv",
                  "SizeDistributionParameters," +
                      std::to_string(FracSizeDistriPara.x) + "," +
                      std::to_string(FracSizeDistriPara.y) + "," +
                      std::to_string(FracSizeDistriPara.z) + "," +
                      std::to_string(FracSizeDistriPara.w) + ",\n");
    AddLineToFile("" + DFNFileName + "/" + DFN_Gen_csv_name +
                      ".csv",
                  "Beta," + std::to_string(Beta) + ",\n");
    AddLineToFile(
        "" + DFNFileName + "/" + DFN_Gen_csv_name + ".csv",
        "Gamma," + DoubleNumberToScientificNotationString(Gamma) +
            ",\n");
};