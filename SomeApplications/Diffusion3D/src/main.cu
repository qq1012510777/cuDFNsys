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
bool IfAFileExist(const string &FilePath);

int main(int argc, char *argv[])
{
    if (argc != 5)
    {
        std::cout << "Usage: " << argv[0] << 
        " <csvDFN> <csvMesh> <csvFlow> <csvPT>\n";
        exit(0);
    }
    std::string csvDFN = std::string(argv[1]);
    std::string csvMesh = std::string(argv[2]);
    std::string csvFlow = std::string(argv[3]);
    std::string csvPT = std::string(argv[4]);

    std::string DFNH5 = "Class_DFN";
    std::string MeshH5 = "Class_Mesh";
    std::string FlowH5 = "Class_Mesh";

    cuDFNsys::DFN<double> my_dfn;
    cuDFNsys::MeshDFN<double> my_mesh;
    cuDFNsys::FlowDFN<double> my_flow;
    cuDFNsys::PTDFN<double> my_PT;

    if (IfAFileExist(DFNH5 + ".h5"))
    {
        my_dfn.LoadDFNFromCSV(csvDFN);
        my_dfn.IdentifyIntersectionsClusters(true);
        std::string DFNVisualFIle = "DFN_Visual";
        my_dfn.Visualization(DFNVisualFIle, DFNVisualFIle, DFNVisualFIle, false, false, true, true);
        my_dfn.StoreInH5(DFNH5);

        my_mesh.LoadParametersFromCSV(csvMesh);
        my_mesh.MeshGeneration(my_dfn);
        std::string MeshVisualFIle = "Mesh_Visual";
        my_mesh.Visualization(my_dfn, MeshVisualFIle, MeshVisualFIle, MeshVisualFIle, false, false);
        my_mesh.StoreInH5(MeshH5);

        my_flow.LoadParametersFromCSV(csvFlow);
        my_flow.FlowSimulation(my_dfn, my_mesh);
        std::string FlowVisualFIle = "Flow_Visual";
        my_flow.Visualization(my_dfn,my_mesh, FlowVisualFIle, FlowVisualFIle, FlowVisualFIle);
        my_flow.FlowData.PressureInteriorEdge = my_flow.FlowData.PressureInteriorEdge.setZero();
        my_flow.FlowData.PressureEles = my_flow.FlowData.PressureEles.setZero();
        my_flow.FlowData.VelocityNormalScalarSepEdges = my_flow.FlowData.VelocityNormalScalarSepEdges.setZero();
        // diffusion no velocity
        my_flow.StoreInH5(FlowH5);
        return 0;
    }
    else
    {
        my_dfn.LoadClassFromH5(DFNH5);
        my_mesh.LoadClassFromH5(MeshH5);
        my_flow.LoadClassFromH5(FlowH5);   
    };


    my_PT.LoadParametersFromCSV(csvPT);
    my_PT.ParticleTracking(my_dfn,
        my_mesh,
        my_flow);
    std::string PTVisualFIle = "PT_Visual";
    my_PT.Visualization(my_dfn,
        my_mesh,
        my_flow, PTVisualFIle, PTVisualFIle, PTVisualFIle);

    return 0;
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