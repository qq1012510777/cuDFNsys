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
bool AreVectorsEqual(const std::vector<size_t> &vec1, const std::vector<size_t> &vec2);

void GetFlowData(const cuDFNsys::DFN<double> &my_dfn,
                 const cuDFNsys::MeshDFN<double> &my_mesh,
                 const cuDFNsys::FlowDFN<double> &my_flow,
                 double &Permeability_apparent,
                 double &q_x,
                 double &q_y,
                 double &q_z);
bool DoesDatasetExist(const std::string &fileName, const std::string &datasetName);

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

    int endII = 0;
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
        if (i == NumGroups - 1)
            endII = 19 + index;
    }

    double minGridSize = atof(argv[endII + 1]);
    endII++;
    double maxGridSize = atof(argv[endII + 1]);
    endII++;

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
    cout << "Mesh size: " << minGridSize << ", " << maxGridSize << endl;
    //--------------------------------------
    int errorCountLimit = 5;
    int erroCounter = 0;
    double P33_total_A[1] = {0};
    double P33_connected_A[1] = {0};
    double Ratio_of_P33_A[1] = {0};
    double P33_largest_cluster_A[1] = {0};
    double P32_total_A[1] = {0};
    double P32_connected_A[1] = {0};
    double Ratio_of_P32_A[1] = {0};
    double P32_largest_cluster_A[1] = {0};
    double P30_A[1] = {0};
    double P30_connected_A[1] = {0};
    double Ratio_of_P30_A[1] = {0};
    double P30_largest_cluster_A[1] = {0};
    double Percolation_probability_A[1] = {0};
    double n_I_A[1] = {0};

    for (int i = 0; i <= NumFracIncre_input[0]; ++i)
    {
        try
        {
            cout << "\n\n---------------loop " << i + 1 << endl;
            /* code */
            cuDFNsys::DFN<double> my_dfn;
            cuDFNsys::MeshDFN<double> my_mesh;

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
                string dataFile = "DataDFN_" + cuDFNsys::ToStringWithWidth(i, 2) + "_" +
                                  cuDFNsys::ToStringWithWidth(percoDir, 1) + ".h5";
                if (IfAFileExist(dataFile))
                {

                    if (DoesDatasetExist(dataFile, "q_z"))
                    {
                        cout << (dataFile)
                             << " exists\n";
                        continue;
                    }
                }
                IfReRun = true;
            }

            // string DFNh5 = "ClassDFN_" + cuDFNsys::ToStringWithWidth(i, 2);
            if (IfReRun)
            {
                cout << "This loop is running\n";

                std::vector<std::pair<bool, std::vector<size_t>>> PercoIndicator(3);

                for (int percoDir = 0; percoDir < 3; ++percoDir)
                {

                    my_dfn.PercoDir = percoDir;

                    if (percoDir == 0)
                    {
                        my_dfn.FractureGeneration();
                        my_dfn.IdentifyIntersectionsClusters(true);
                        // my_dfn.StoreInH5(DFNh5);
                    }
                    else
                    {
                        my_dfn.ChangePercolationDirectionIdentifyPercolationCluster(percoDir);
                    }
                    std::cout << "PercolationDirection = " << my_dfn.PercoDir << ", PercolationClusterSize = " << my_dfn.PercolationCluster.size() << endl;
                    PercoIndicator[percoDir].first = (my_dfn.PercolationCluster.size() > 0 ? true : false);
                    if (PercoIndicator[percoDir].first)
                        PercoIndicator[percoDir].second = my_dfn.PercolationCluster;

                    cuDFNsys::GetStatistics<double>(my_dfn.FracturesHost,
                                                    my_dfn.IntersectionMap,
                                                    my_dfn.ListClusters,
                                                    my_dfn.PercolationCluster,
                                                    my_dfn.DomainSizeX,
                                                    P33_total_A[0],
                                                    P33_connected_A[0],
                                                    Ratio_of_P33_A[0],
                                                    P33_largest_cluster_A[0],
                                                    P32_total_A[0],
                                                    P32_connected_A[0],
                                                    Ratio_of_P32_A[0],
                                                    P32_largest_cluster_A[0],
                                                    P30_A[0],
                                                    P30_connected_A[0],
                                                    Ratio_of_P30_A[0],
                                                    P30_largest_cluster_A[0],
                                                    Percolation_probability_A[0],
                                                    n_I_A[0]);

                    vector<string> datasetname = {
                        "P33_total",
                        "P33_connected",
                        "Ratio_of_P33",
                        "P33_largest_cluster",
                        "P32_total",
                        "P32_connected",
                        "Ratio_of_P32",
                        "P32_largest_cluster",
                        "P30",
                        "P30_connected",
                        "Ratio_of_P30",
                        "P30_largest_cluster",
                        "Percolation_probability",
                        "n_I"};

                    vector<double *> data_input = {P33_total_A,
                                                   P33_connected_A,
                                                   Ratio_of_P33_A,
                                                   P33_largest_cluster_A,
                                                   P32_total_A,
                                                   P32_connected_A,
                                                   Ratio_of_P32_A,
                                                   P32_largest_cluster_A,
                                                   P30_A,
                                                   P30_connected_A,
                                                   Ratio_of_P30_A,
                                                   P30_largest_cluster_A,
                                                   Percolation_probability_A,
                                                   n_I_A

                    };
                    string dataFile = "DataDFN_" + cuDFNsys::ToStringWithWidth(i, 2) + "_" +
                                      cuDFNsys::ToStringWithWidth(percoDir, 1) + ".h5";

                    vector<uint2> dim_ss(data_input.size(), make_uint2(1, 1));
                    cuDFNsys::HDF5API h5out;
                    h5out.NewFile(dataFile);
                    h5out.AddDatasetsWithOneGroup(dataFile, "N",
                                                  datasetname, data_input, dim_ss);
                }

                bool diffMesh = false;
                for (int j = 0; j < 3; ++j)
                    for (int k = 0; k < 3; ++k)
                        if (!AreVectorsEqual(PercoIndicator[j].second, PercoIndicator[k].second))
                            diffMesh = true;

                cuDFNsys::FlowDFN<double> my_flow;
                cout << "\n";
                cout << "If generate different mesh: " << (diffMesh ? "yes" : "no") << endl;
                cout << "\n";
                for (int percoDir = 0; percoDir < 3; ++percoDir)
                {
                    // if (i == 4)
                    //     cout << "Warning\n";
                    //
                    string dataFile = "DataDFN_" + cuDFNsys::ToStringWithWidth(i, 2) + "_" +
                                      cuDFNsys::ToStringWithWidth(percoDir, 1) + ".h5";

                    double q_x = 0, q_y = 0, q_z = 0;
                    double Permeability_apparent = 0;

                    std::cout << "PercolationDirection = " << percoDir << ", flowSimu ...\n";

                    if (PercoIndicator[percoDir].first)
                    {
                        if (!diffMesh) // do not need generate new mesh
                        {
                            // for dense networks, all directions are percolative directions
                            my_dfn.ChangePercolationDirectionIdentifyPercolationCluster(percoDir);
                            if (percoDir == 0)
                            {
                                my_mesh.MinElementSize = minGridSize;
                                my_mesh.MaxElementSize = maxGridSize;
                                my_mesh.MeshGeneration(my_dfn);
                            }
                            else
                            {
                                my_mesh.ChangePecolationDirectionAndRenumberingEdge(percoDir,
                                                                                    my_dfn.DomainSizeX,
                                                                                    my_dfn.DomainDimensionRatio);
                            }
                            //----------simuation flow
                            my_flow.InletHead = my_dfn.DomainSizeX;
                            my_flow.OutletHead = 0;

                            my_flow.FlowSimulation(my_dfn, my_mesh);
                            GetFlowData(my_dfn,
                                        my_mesh,
                                        my_flow,
                                        Permeability_apparent,
                                        q_x,
                                        q_y,
                                        q_z);
                            // my_flow.StoreInH5("DFNFlow_" + cuDFNsys::ToStringWithWidth(i, 2) + "_" +
                            //                   cuDFNsys::ToStringWithWidth(percoDir, 1) + ".h5");
                            // my_flow.Visualization(my_dfn, my_mesh,
                            //                       "MatlabDFNFlow_" + cuDFNsys::ToStringWithWidth(i, 2) + "_" +
                            //                           cuDFNsys::ToStringWithWidth(percoDir, 1),
                            //                       "MatlabDFNFlow_" + cuDFNsys::ToStringWithWidth(i, 2) + "_" +
                            //                           cuDFNsys::ToStringWithWidth(percoDir, 1),
                            //                       "MatlabDFNFlow_" + cuDFNsys::ToStringWithWidth(i, 2) + "_" +
                            //                           cuDFNsys::ToStringWithWidth(percoDir, 1));
                        }
                        else
                        {
                            // for sparse networks, maybe percolation cluster is different
                            cuDFNsys::DFN<double> my_dfn_2 = my_dfn;
                            my_dfn_2.ChangePercolationDirectionIdentifyPercolationCluster(percoDir);

                            cuDFNsys::MeshDFN<double> my_mesh_3;
                            my_mesh_3.MinElementSize = minGridSize;
                            my_mesh_3.MaxElementSize = maxGridSize;
                            my_mesh_3.MeshGeneration(my_dfn_2);

                            //----------simuation flow
                            my_flow.InletHead = my_dfn_2.DomainSizeX;
                            my_flow.OutletHead = 0;
                            my_flow.FlowSimulation(my_dfn_2, my_mesh_3);

                            GetFlowData(my_dfn_2,
                                        my_mesh_3,
                                        my_flow,
                                        Permeability_apparent,
                                        q_x,
                                        q_y,
                                        q_z);
                            // my_flow.StoreInH5("DFNFlow_" + cuDFNsys::ToStringWithWidth(i, 2) + "_" +
                            //                   cuDFNsys::ToStringWithWidth(percoDir, 1) + ".h5");
                            // my_flow.Visualization(my_dfn_2, my_mesh_3,
                            //                       "MatlabDFNFlow_" + cuDFNsys::ToStringWithWidth(i, 2) + "_" +
                            //                           cuDFNsys::ToStringWithWidth(percoDir, 1),
                            //                       "MatlabDFNFlow_" + cuDFNsys::ToStringWithWidth(i, 2) + "_" +
                            //                           cuDFNsys::ToStringWithWidth(percoDir, 1),
                            //                       "MatlabDFNFlow_" + cuDFNsys::ToStringWithWidth(i, 2) + "_" +
                            //                           cuDFNsys::ToStringWithWidth(percoDir, 1));
                        }
                    }

                    cuDFNsys::HDF5API h5out;
                    h5out.AddDataset(dataFile, "N", "PermeabilityApparent", &Permeability_apparent, make_uint2(1, 1));
                    h5out.AddDataset(dataFile, "N", "q_x", &q_x, make_uint2(1, 1));
                    h5out.AddDataset(dataFile, "N", "q_y", &q_y, make_uint2(1, 1));
                    h5out.AddDataset(dataFile, "N", "q_z", &q_z, make_uint2(1, 1));

                    // if (!diffMesh && PercoIndicator[percoDir].first)
                    //     exit(0);
                }
            }

            erroCounter = 0;
        }
        catch (cuDFNsys::ExceptionsIgnore e)
        {
            cout << "cuDFNsys::Exceptions: " << e.what() << endl;
            int result_system = system(
                std::string("rm -rf ./ClassDFN_" + cuDFNsys::ToStringWithWidth(i, 2) + ".h5")
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
                std::string("rm -rf ./ClassDFN_" + cuDFNsys::ToStringWithWidth(i, 2) + ".h5")
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
                std::string("rm -rf ./ClassDFN_" + cuDFNsys::ToStringWithWidth(i, 2) + ".h5")
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
                std::string("rm -rf ./ClassDFN_" + cuDFNsys::ToStringWithWidth(i, 2) + ".h5")
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
                std::string("rm -rf ./ClassDFN_" + cuDFNsys::ToStringWithWidth(i, 2) + ".h5")
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

bool AreVectorsEqual(const std::vector<size_t> &vec1, const std::vector<size_t> &vec2)
{
    // Check if sizes are equal
    if (vec1.size() != vec2.size())
    {
        return false;
    }

    // Compare element by element
    return std::equal(vec1.begin(), vec1.end(), vec2.begin());
}
void GetFlowData(const cuDFNsys::DFN<double> &my_dfn,
                 const cuDFNsys::MeshDFN<double> &my_mesh,
                 const cuDFNsys::FlowDFN<double> &my_flow,
                 double &Permeability_apparent,
                 double &q_x,
                 double &q_y,
                 double &q_z)
{
    //-----------------get data
    Permeability_apparent = my_flow.FlowData.Permeability;

    for (uint j = 0; j < my_mesh.MeshData.Element3D.size(); ++j)
    {
        uint3 EdgeNO = make_uint3((j + 1) * 3 - 3, (j + 1) * 3 - 2,
                                  (j + 1) * 3 - 1); // from 0
        cuDFNsys::Vector3<double> Velocity_ = cuDFNsys ::MakeVector3(
            my_flow.FlowData.VelocityNormalScalarSepEdges(EdgeNO.x, 0),
            my_flow.FlowData.VelocityNormalScalarSepEdges(EdgeNO.y, 0),
            my_flow.FlowData.VelocityNormalScalarSepEdges(EdgeNO.z, 0));
        cuDFNsys::Vector2<double> Vertexes[3];
        Vertexes[0] = cuDFNsys::MakeVector2(my_mesh.MeshData.Coordinate2D[j].x[0],
                                            my_mesh.MeshData.Coordinate2D[j].y[0]);
        Vertexes[1] = cuDFNsys::MakeVector2(my_mesh.MeshData.Coordinate2D[j].x[1],
                                            my_mesh.MeshData.Coordinate2D[j].y[1]);
        Vertexes[2] = cuDFNsys::MakeVector2(my_mesh.MeshData.Coordinate2D[j].x[2],
                                            my_mesh.MeshData.Coordinate2D[j].y[2]);

        cuDFNsys::Vector2<double> Center_p = cuDFNsys::MakeVector2(
            1.0f / 3.0f * (Vertexes[0].x + Vertexes[1].x + Vertexes[2].x),
            1.0f / 3.0f * (Vertexes[0].y + Vertexes[1].y + Vertexes[2].y));

        cuDFNsys::Vector2<double> velocity_p =
            cuDFNsys::ReconstructVelocityGrid<double>(Center_p, Vertexes, Velocity_);

        cuDFNsys::Vector3<double> velocity_p_3D =
            cuDFNsys::MakeVector3(velocity_p.x, velocity_p.y, 0.);

        double R_mat[3][3];
        cuDFNsys::Fracture<double> FII =
            my_dfn.FracturesHost[my_mesh.MeshData.ElementFracTag[j]];
        FII.RoationMatrix(R_mat, 23);
        velocity_p_3D =
            cuDFNsys::ProductSquare3Vector3<double>(R_mat, velocity_p_3D);

        // double b_f = pow(FII.Conductivity * 12, 1.0 / 3.0);
        // velocity_p_3D.x /= b_f,
        //     velocity_p_3D.y /= b_f,
        //     velocity_p_3D.z /= b_f;

        double area_this_element =
            cuDFNsys::Triangle2DArea<double>(Vertexes[0],
                                             Vertexes[1],
                                             Vertexes[2]);
        // cout << j + 1 << ":  " << velocity_p_3D << ", a: " << area_this_element << endl;
        velocity_p_3D.x *= area_this_element,
            velocity_p_3D.y *= area_this_element,
            velocity_p_3D.z *= area_this_element;
        q_x += velocity_p_3D.x,
            q_y += velocity_p_3D.y,
            q_z += velocity_p_3D.z;
    }

    q_x /= pow(my_dfn.DomainSizeX, 3.0);
    q_y /= pow(my_dfn.DomainSizeX, 3.0);
    q_z /= pow(my_dfn.DomainSizeX, 3.0);
}
bool DoesDatasetExist(const std::string &fileName, const std::string &datasetName)
{
    try
    {
        // Open the file in read-only mode
        H5::H5File file(fileName, H5F_ACC_RDONLY);

        // Check if the dataset exists
        if (H5Lexists(file.getId(), datasetName.c_str(), H5P_DEFAULT) > 0)
        {
            return true;
        }
    }
    catch (const H5::Exception &e)
    {
        std::cerr << "HDF5 Error: " << e.getDetailMsg() << std::endl;
    }
    return false;
}