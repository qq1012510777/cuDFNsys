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

std::vector<double> HeadGradient = {1};
std::vector<double> Lm = {DomainSizeX};
std::vector<double> Num_Fracs = {0};
std::vector<double> P32_LargestCluster = {0};
std::vector<double> P32_Total = {0};
std::vector<double> P33_LargestCluster = {0};
std::vector<double> P33_Total = {0};
std::vector<double> Percolation_Status = {0};
std::vector<double> Permeability_Apparent = {0};
std::vector<double> q = {0, 0, 0};

string GenCsvFiles(const string &DFNFileName, const int &stepNo_i);
void CreateOrEmptyFile(const std::string &filename);
void AddLineToFile(const std::string &filename, const std::string &line);
bool IfAFileExist(const string &FilePath);
string DoubleNumberToScientificNotationString(const double &number);
void RecordDataH5(const string &FileNameH5RecordData);

int main(int argc, char *argv[])
{

    string cwd = fs::current_path().string();

    time_t t;
    time(&t);

    ExeuctablePath = argv[1];
    // DFN generation parameters
    DomainSizeX = atof(argv[2]);
    DomainDimensionRatio =
        make_double3(atof(argv[3]), atof(argv[4]), atof(argv[5]));
    // PercoDir = 0;
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
        // string DFNFileName = cwd + "/DFN_" + cuDFNsys::ToStringWithWidth(i + 1, 3);
        // int result_system = system(("mkdir " + DFNFileName + " -p").c_str());

        // --------------- first we generate a DFN
        string csv_DFN_gen = GenCsvFiles(cwd, i);

        cuDFNsys::DFN<double> dfnGen_ori;
        dfnGen_ori.LoadDFNFromCSV(csv_DFN_gen);
        string original_DFN_H5 = "Class_DFN_orignal_" + cuDFNsys::ToStringWithWidth(i + 1, 3);
        dfnGen_ori.StoreInH5(original_DFN_H5);

        //-----------------
        // some data should be recorded
        HeadGradient = {1};
        Lm = {DomainSizeX};
        Num_Fracs = {0};
        P32_LargestCluster = {0};
        P32_Total = {0};
        P33_LargestCluster = {0};
        P33_Total = {0};
        Percolation_Status = {0};
        Permeability_Apparent = {0};
        q = {0, 0, 0};

        for (int j = 0; j < 3; ++j)
        {
            string DFNdataFileNameH5 = "data_" + cuDFNsys::ToStringWithWidth(i + 1, 3) + ".h5";
            //-----------------------
            cuDFNsys::DFN<double> dfnGen;

            try
            {
                dfnGen.LoadClassFromH5(original_DFN_H5);
                dfnGen.PercoDir = j;
                dfnGen.IdentifyIntersectionsClusters(false);

                Num_Fracs[0] = dfnGen.FracturesHost.size();

                for (const auto &e : dfnGen.FracturesHost)
                    P32_Total[0] += pow(sqrt(2.) * (e.Radius), 2.),
                        P33_Total[0] += pow(sqrt(2.) * (e.Radius), 2.) * pow(e.Conductivity * 12., 1. / 3.);
                P32_Total[0] /= pow(DomainSizeX, 3.),
                    P33_Total[0] /= pow(DomainSizeX, 3.);

                for (size_t k = 0; k < dfnGen.PercolationCluster.size(); ++k)
                    for (size_t l = 0; l < dfnGen.ListClusters[dfnGen.PercolationCluster[k]].size(); ++l)
                        P32_LargestCluster[0] += pow(sqrt(2.) * (dfnGen.FracturesHost[dfnGen.ListClusters[dfnGen.PercolationCluster[k]][l]].Radius), 2.),
                            P33_LargestCluster[0] += pow(sqrt(2.) * (dfnGen.FracturesHost[dfnGen.ListClusters[dfnGen.PercolationCluster[k]][l]].Radius), 2.) *
                                                     pow(dfnGen.FracturesHost[dfnGen.ListClusters[dfnGen.PercolationCluster[k]][l]].Conductivity * 12., 1. / 3.);

                dfnGen.IdentifyIntersectionsClusters(true);
            }
            catch (...)
            {
            }

            if (dfnGen.PercolationCluster.size() == 0)
                RecordDataH5(DFNdataFileNameH5);
            continue;
        }
    }

    return 0;
}

string GenCsvFiles(const string &DFNFileName, const int &stepNo_i)
{
    string DFN_Gen_csv_name = "StochasticDFN_" + cuDFNsys::ToStringWithWidth(stepNo_i + 1, 3);

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
    return DFNFileName + "/" + DFN_Gen_csv_name;
};
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
string DoubleNumberToScientificNotationString(const double &number)
{
    std::stringstream ss;
    ss << std::scientific << number;
    std::string result = ss.str();
    return result;
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
void RecordDataH5(const string &FileNameH5RecordData)
{
    cuDFNsys::HDF5API h5g;
    h5g.NewFile(FileNameH5RecordData);
    h5g.AddDataset(FileNameH5RecordData, "/", "HeadGradient", tmHeadGradientpA.data(), make_uint2{1, 1});

    h5g.AddDataset(FileNameH5RecordData, "/", "Lm", Lm.data(), make_uint2{1, 1});
    h5g.AddDataset(FileNameH5RecordData, "/", "Num_Fracs", Num_Fracs.data(), make_uint2{1, 1});
    h5g.AddDataset(FileNameH5RecordData, "/", "P32_LargestCluster", P32_LargestCluster.data(), make_uint2{1, 1});
    h5g.AddDataset(FileNameH5RecordData, "/", "P32_Total", P32_Total.data(), make_uint2{1, 1});
    h5g.AddDataset(FileNameH5RecordData, "/", "P33_LargestCluster", P33_LargestCluster.data(), make_uint2{1, 1});
    h5g.AddDataset(FileNameH5RecordData, "/", "P33_Total", P33_Total.data(), make_uint2{1, 1});
    h5g.AddDataset(FileNameH5RecordData, "/", "Percolation_Status", Percolation_Status.data(), make_uint2{1, 1});
    h5g.AddDataset(FileNameH5RecordData, "/", "Permeability_Apparent", Permeability_Apparent.data(), make_uint2{1, 1});
    h5g.AddDataset(FileNameH5RecordData, "/", "q", q.data(), make_uint2{3, 1});
};