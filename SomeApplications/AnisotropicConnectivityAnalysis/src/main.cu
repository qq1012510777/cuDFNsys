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

std::ostream &operator<<(std::ostream &os, const double3 &vec)
{
    os << vec.x << ", " << vec.y << ", " << vec.z;
    return os;
}
std::ostream &operator<<(std::ostream &os, const double4 &vec)
{
    os << vec.x << ", " << vec.y << ", " << vec.z << ", " << vec.w;
    return os;
}

bool ErrorCheck(int &erroCounter, const int &errorCountLimit)
{
    erroCounter++;
    if (erroCounter > errorCountLimit)
    {
        std::cout << "Too many errors happen\n";
        return true;
    }
    return false;
};

int main(int argc, char *argv[])
{
    time_t t;
    time(&t);

    // string ExeuctablePath = argv[1];
    double DomainSizeX = atof(argv[1]);
    double3 DomainDimensionRatio =
        make_double3(atof(argv[2]), atof(argv[3]), atof(argv[4]));

    int NumGroups = atoi(argv[5]);

    std::vector<double> Kappa_input(NumGroups, 0);
    std::vector<cuDFNsys::Vector3<double>> MeanOrientationOfFisherDistribution_input(NumGroups);
    std::vector<int> FracNumInit_input(NumGroups, 0),
        FracNumIncre_input(NumGroups, 0),
        NumFracIncre_input(NumGroups, 0);
    std::vector<int> ModeOfSizeDistribution_input(NumGroups, 0);
    std::vector<cuDFNsys::Vector4<double>> SizeDistributionParameters_input(NumGroups);
    std::vector<double> Beta_input(NumGroups, 0), Gamma_input(NumGroups, 0);

    for (int i = 0; i < NumGroups; ++i)
    {
        int index = 14 * i;

        Kappa_input[i] = atof(argv[6 + index]);
        MeanOrientationOfFisherDistribution_input[i] = make_double3(
            atof(argv[7 + index]), atof(argv[8 + index]), atof(argv[9 + index]));
        FracNumInit_input[i] = atoi(argv[10 + index]);
        FracNumIncre_input[i] = atoi(argv[11 + index]);
        NumFracIncre_input[i] = atoi(argv[12 + index]);

        ModeOfSizeDistribution_input[i] = atoi(argv[13 + index]);

        SizeDistributionParameters_input[i] = make_double4(
            atof(argv[14 + index]), atof(argv[15 + index]),
            atof(argv[16 + index]), atof(argv[17 + index]));

        Beta_input[i] = atof(argv[18 + index]);
        Gamma_input[i] = atof(argv[19 + index]);
    }
    std::cout << "Domain Size = " << DomainSizeX << endl;
    cout << "Domain dimension ratio: " << DomainDimensionRatio << endl;
    cout << "Number of fracture groups: " << NumGroups << endl;

    for (int i = 0; i < NumGroups; ++i)
    {
        cout << "---------group " << i + 1 << "----------------\n";
        cout << "Kappa: " << Kappa_input[i] << endl;
        cout << "Mean orientation: " << MeanOrientationOfFisherDistribution_input[i] << endl;
        cout << "Init fracture number: " << FracNumInit_input[i] << endl;
        cout << "Facture number increament: " << FracNumIncre_input[i] << endl;
        cout << "Number of increments: " << NumFracIncre_input[i] << endl;
        cout << "ModeOfSizeDistribution_input: " << ModeOfSizeDistribution_input[i] << endl;
        cout << "SizeDistributionParameters_input: " << SizeDistributionParameters_input[i] << endl;
        cout << "Beta_input: " << Beta_input[i] << endl;
        cout << "Gamma_input: " << Gamma_input[i] << endl;
        cout << "--------------------------------\n\n";
    }
    //--------------------------------------
    int errorCountLimit = 5;
    int erroCounter = 0;

    for (int i = 0; i <= NumFracIncre_input[0]; ++i)
    {
        try
        {
            cout << "\n\n---------------loop " << i + 1 << endl;
            /* code */
            cuDFNsys::DFN<double> my_dfn;
            my_dfn.RandomSeed = (unsigned long)t;
            my_dfn.DomainSizeX = DomainSizeX;
            my_dfn.DomainDimensionRatio = DomainDimensionRatio;

            my_dfn.NumFractures.resize(NumGroups);
            my_dfn.Kappa.resize(NumGroups);
            my_dfn.MeanOrientationOfFisherDistribution.resize(NumGroups);
            my_dfn.Beta.resize(NumGroups);
            my_dfn.Gamma.resize(NumGroups);
            my_dfn.ModeOfSizeDistribution.resize(NumGroups);
            my_dfn.SizeDistributionParameters.resize(NumGroups);

            for (int j = 0; j < NumGroups; ++j)
            {
                my_dfn.NumFractures[j] = FracNumInit_input[j] + i * FracNumIncre_input[j];
                my_dfn.Kappa[j] = Kappa_input[j];
                my_dfn.MeanOrientationOfFisherDistribution[j] = MeanOrientationOfFisherDistribution_input[j];
                my_dfn.Beta[j] = Beta_input[j];
                my_dfn.Gamma[j] = Gamma_input[j];
                my_dfn.ModeOfSizeDistribution[j] = ModeOfSizeDistribution_input[j];
                my_dfn.SizeDistributionParameters[j] = SizeDistributionParameters_input[j];
            }

            bool IfReRun = false;
            for (int percoDir = 0; percoDir < 3; ++percoDir)
            {
                if (IfAFileExist("ClassDFN_" + cuDFNsys::ToStringWithWidth(i, 2) + "_" +
                                 cuDFNsys::ToStringWithWidth(percoDir, 1) + ".h5"))
                {
                    cout << ("ClassDFN_" + cuDFNsys::ToStringWithWidth(i, 2) + "_" +
                             cuDFNsys::ToStringWithWidth(percoDir, 1) + ".h5")
                         << " exists\n";
                    continue;
                }
                IfReRun = true;
            }

            if (IfReRun)
            {
                cout << "This loop is running\n";
                for (int percoDir = 0; percoDir < 3; ++percoDir)
                {

                    cuDFNsys::DFN<double> my_dfn_2;
                    my_dfn_2 = my_dfn;

                    my_dfn_2.PercoDir = percoDir;
                    my_dfn_2.FractureGeneration();
                    my_dfn_2.IdentifyIntersectionsClusters(false);

                    my_dfn_2.StoreInH5("ClassDFN_" + cuDFNsys::ToStringWithWidth(i, 2) + "_" +
                                       cuDFNsys::ToStringWithWidth(percoDir, 1));
                }
            }

            erroCounter = 0;
        }
        catch (cuDFNsys::ExceptionsIgnore e)
        {
            cout << "cuDFNsys::Exceptions: " << e.what() << endl;
            int result_system = system(
                std::string("rm -rf ./ClassDFN_" + cuDFNsys::ToStringWithWidth(i, 2) + "_*.h5")
                    .c_str());

            i--;
            if (ErrorCheck(erroCounter, errorCountLimit))
                break;
            continue;
        }
        catch (cuDFNsys::ExceptionsPause e)
        {
            cout << "cuDFNsys::Exceptions: " << e.what() << endl;
            int result_system = system(
                std::string("rm -rf ./ClassDFN_" + cuDFNsys::ToStringWithWidth(i, 2) + "_*.h5")
                    .c_str());
            i--;
            if (ErrorCheck(erroCounter, errorCountLimit))
                break;
            continue;
        }
        catch (H5::Exception e)
        {
            cout << "H5::Exceptions: " << e.getDetailMsg() << endl;
            int result_system = system(
                std::string("rm -rf ./ClassDFN_" + cuDFNsys::ToStringWithWidth(i, 2) + "_*.h5")
                    .c_str());
            i--;
            if (ErrorCheck(erroCounter, errorCountLimit))
                break;
            continue;
        }
        catch (H5::FileIException e)
        {
            cout << "H5::Exceptions: " << e.getDetailMsg() << endl;
            int result_system = system(
                std::string("rm -rf ./ClassDFN_" + cuDFNsys::ToStringWithWidth(i, 2) + "_*.h5")
                    .c_str());
            i--;
            if (ErrorCheck(erroCounter, errorCountLimit))
                break;
            continue;
        }
        catch (...)
        {
            cout << "Unknown exceptions\n";
            int result_system = system(
                std::string("rm -rf ./ClassDFN_" + cuDFNsys::ToStringWithWidth(i, 2) + "_*.h5")
                    .c_str());
            i--;
            if (ErrorCheck(erroCounter, errorCountLimit))
                break;
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
