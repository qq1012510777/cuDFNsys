
#include "cuDFNsys.cuh"
#include <filesystem>
#include <fstream>

namespace fs = std::filesystem;

const int NumDataPnts = 30;
// we consider a uniformly-orientated DFN with mono-sized fractures
// so the radius of fractures (r) is set to be 7.5 m
// Fisher constant kappa = 0
// dimensionless density = density of fractures times a local volume scale
// the value of the local volume scale is `V_ex_a` here
// dimensionless density ranges from 1 to 11.7300
// the maximum densionless density is around 5 times the percolation threshold (2.31)
const double V_ex_a = 2.386485386504598e+03;
const double Chosen_rho_prime_value[NumDataPnts] = {
    1.0000, 1.3700,  1.7400,  2.1100,  2.4800,  2.8500, 3.2200, 3.5900,
    3.9600, 4.3300,  4.7000,  5.0700,  5.4400,  5.8100, 6.1800, 6.5500,
    6.9200, 7.2900,  7.6600,  8.0300,  8.4000,  8.7700, 9.1400, 9.5100,
    9.8800, 10.2500, 10.6200, 10.9900, 11.3600, 11.7300};
const double AreaOfAFracture = 112.5000;

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

    my_dfn.RandomSeed = (unsigned long)t * i + MC_NO_;
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
        meshGen.MinElementSize = 0.5;
        meshGen.MaxElementSize = 1.;

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

    h5g.AddDataset(DataFileName, GroupName, "Conn", &Conn, make_uint2(1, 0));
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
    return 0;
};
