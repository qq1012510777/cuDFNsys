
#include "cuDFNsys.cuh"
#include <H5Exception.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

const int NumDataPnts = 30;
// we consider a uniformly-orientated DFN with mono-sized fractures
// so the radius of fractures (r) is set to be 7.5 m
// Fisher constant kappa = 0
// dimensionless density = density of fractures times a local volume scale
// the value of the local volume scale is `V_ex_a` here
const double V_ex_a = 2.386485386504598e+03;
const double AreaOfAFracture = 112.5000;

std::vector<double> readCSV(const std::string &filename);

int main(int argc, char *argv[])
{
    //-----------if a HDF5 file exists? --------
    std::string filename = "Data.h5";
    if (fs::exists(filename))
    {
        std::cout << "File exists in the current directory." << std::endl;

        // check if this has data

        cuDFNsys::HDF5API TRES;
        bool If_has_data = true;
        try
        {
            std::vector<double> TMP_v =
                TRES.ReadDataset<double>("Data.h5", "N", "NumFractures");
        }
        catch (...)
        {
            If_has_data = false;
        }

        if (If_has_data)
        {
            std::ofstream outputFile;
            outputFile.open("Finished");
            outputFile.close();
            return 0;
        }
    }
    else
    {
        std::cout << "File does not exist in the current directory."
                  << std::endl;
    }
    //--------------------

    //----------read rho^prime data

    std::vector<double> Chosen_rho_prime_value =
        readCSV("../../Chosen_rho_prime.csv");
    if (Chosen_rho_prime_value.size() != NumDataPnts)
    {
        cout << "Do not have enough data\n";
        exit(0);
    }

    time_t t;
    time(&t);

    double DomainSize = atof(argv[1]);
    int i = atoi(argv[2]);
    int MC_NO_ = atoi(argv[3]);

    cuDFNsys::HDF5API h5g;
    string DataFileName = "Data.h5";
    h5g.NewFile(DataFileName);

    double Conn = 0, Permeab = 0, NumFiniteElements = 0, DFNTime = 0,
           MeshTime = 0, FlowTime = 0, NumFractures = 0;

    int NumFracs_i =
        round(Chosen_rho_prime_value[i] / V_ex_a * pow(DomainSize, 3));

    NumFractures = NumFracs_i;

    cuDFNsys::DFN<double> my_dfn;
    cuDFNsys::MeshDFN<double> meshGen;
    cuDFNsys::FlowDFN<double> flowDFN;

    try
    {
        my_dfn.RandomSeed = (unsigned long)(t + i + MC_NO_);

        {
            // record random seed for reproducing, if error happens
            h5g.NewFile("RandomSeed.h5");
            double ls = my_dfn.RandomSeed;
            h5g.AddDataset("RandomSeed.h5", "N", "RandomSeed", &ls,
                           make_uint2(1, 0));
        }

        my_dfn.NumFractures = {NumFracs_i};
        my_dfn.Kappa = {0};
        my_dfn.MeanOrientationOfFisherDistribution = {make_double3(0, 0, 1)};
        my_dfn.DomainSizeX = DomainSize;
        my_dfn.DomainDimensionRatio = {make_double3(1, 1, 1)};
        my_dfn.Beta = {0.2};
        my_dfn.Gamma = {5e-6};
        my_dfn.ModeOfSizeDistribution = {3};
        my_dfn.SizeDistributionParameters = {make_double4(7.5, 0, 0, 0)};
        my_dfn.PercoDir = 2;

        double i_start_time = cuDFNsys::CPUSecond();

        my_dfn.FractureGeneration();
        my_dfn.IdentifyIntersectionsClusters(false);

        DFNTime = cuDFNsys::CPUSecond() - i_start_time;

        if (my_dfn.PercolationCluster.size() > 0)
        {
            for (auto e : my_dfn.PercolationCluster)
                Conn += (my_dfn.ListClusters[e].size() * AreaOfAFracture);

            Conn = Conn / (NumFracs_i * AreaOfAFracture);
        }

        my_dfn.IdentifyIntersectionsClusters(true);

        if (my_dfn.PercolationCluster.size() > 0)
        {
            meshGen.MinElementSize = 1;
            meshGen.MaxElementSize = 3.;

            i_start_time = cuDFNsys::CPUSecond();

            meshGen.MeshGeneration(my_dfn);

            MeshTime = cuDFNsys::CPUSecond() - i_start_time;

            NumFiniteElements = meshGen.MeshData.Element3D.size();

            flowDFN.MuOverRhoG = 1;
            flowDFN.InletHead = 100;
            flowDFN.OutletHead = 20;

            i_start_time = cuDFNsys::CPUSecond();

            flowDFN.FlowSimulation(my_dfn, meshGen);

            FlowTime = cuDFNsys::CPUSecond() - i_start_time;

            Permeab = flowDFN.FlowData.Permeability;
        }

        string GroupName = "N";

        h5g.AddDataset(DataFileName, GroupName, "Conn", &Conn,
                       make_uint2(1, 0));
        h5g.AddDataset(DataFileName, GroupName, "Permeab", &Permeab,
                       make_uint2(1, 0));
        h5g.AddDataset(DataFileName, GroupName, "NumFiniteElements",
                       &NumFiniteElements, make_uint2(1, 0));

        h5g.AddDataset(DataFileName, GroupName, "DFNTime", &DFNTime,
                       make_uint2(1, 0));
        h5g.AddDataset(DataFileName, GroupName, "MeshTime", &MeshTime,
                       make_uint2(1, 0));
        h5g.AddDataset(DataFileName, GroupName, "FlowTime", &FlowTime,
                       make_uint2(1, 0));

        h5g.AddDataset(DataFileName, GroupName, "NumFractures", &NumFractures,
                       make_uint2(1, 0));

        std::ofstream outputFile;
        outputFile.open("Finished");
        outputFile.close();
    }
    catch (cuDFNsys::ExceptionsIgnore &e)
    {
        cout << e.what() << endl;
        my_dfn.StoreInH5("Class_DFN");
        meshGen.StoreInH5("Class_MESH");
        flowDFN.StoreInH5("Class_FLOW");
    }
    catch (cuDFNsys::ExceptionsPause &e)
    {
        cout << e.what() << endl;
        my_dfn.StoreInH5("Class_DFN");
        meshGen.StoreInH5("Class_MESH");
        flowDFN.StoreInH5("Class_FLOW");
    }
    catch (H5::Exception &e)
    {
        cout << "H5::Exception\n";
        my_dfn.StoreInH5("Class_DFN");
        meshGen.StoreInH5("Class_MESH");
        flowDFN.StoreInH5("Class_FLOW");
    }
    catch (...)
    {
        cout << "Unknown exceptions!\n";
        my_dfn.StoreInH5("Class_DFN");
        meshGen.StoreInH5("Class_MESH");
        flowDFN.StoreInH5("Class_FLOW");
    }
    return 0;
};

std::vector<double> readCSV(const std::string &filename)
{
    std::vector<double> numbers;
    numbers.reserve(NumDataPnts);

    std::ifstream file(filename);
    std::string line;

    if (!file.is_open())
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        return numbers;
    }

    // Read only one line from the file
    if (std::getline(file, line))
    {
        std::stringstream ss(line);
        std::string token;

        while (std::getline(ss, token, ','))
        {
            // Convert the token to a double and add it to the numbers vector
            if (token == "")
                break;
            double value = std::stod(token);
            numbers.push_back(value);
        }
    }

    file.close();
    return numbers;
}
