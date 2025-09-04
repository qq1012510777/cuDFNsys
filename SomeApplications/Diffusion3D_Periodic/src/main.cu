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
    std::string FlowH5 = "Class_Flow";

    cuDFNsys::DFN<double> my_dfn;
    cuDFNsys::MeshDFN<double> my_mesh;
    cuDFNsys::FlowDFN<double> my_flow;
    cuDFNsys::PTDFN<double> my_PT;

    std::string FlowVisualFIle = "Flow_Visual";
    try
    {
        if (!IfAFileExist(DFNH5 + ".h5"))
        {
            my_dfn.LoadDFNFromCSV(csvDFN);
            my_dfn.IdentifyIntersectionsClusters(true);
            if (my_dfn.PercolationCluster.size() == 0)
                return 0;
            std::cout << "**It is percolative**\n";
            std::string DFNVisualFIle = "DFN_Visual";
            my_dfn.SpatialPeriodicity();
            my_dfn.IdentifyIntersectionsClusters(true);

            if (my_dfn.PercolationCluster.size() == 0)
            {
                std::cout << "Wrong with `SpatialPeriodicity`\n";
                return 0;
            }
            my_dfn.Visualization(DFNVisualFIle, DFNVisualFIle, DFNVisualFIle, true, true, true, true);
            my_mesh.LoadParametersFromCSV(csvMesh);
            my_mesh.MeshGeneration(my_dfn);
            std::string MeshVisualFIle = "Mesh_Visual";
            my_mesh.Visualization(my_dfn, MeshVisualFIle, MeshVisualFIle, MeshVisualFIle, true, true);
            
            my_flow.IfPeriodic = true;
            my_flow.LoadParametersFromCSV(csvFlow);
            my_flow.FlowSimulation(my_dfn, my_mesh);
            
            my_flow.Visualization(my_dfn,my_mesh, FlowVisualFIle, FlowVisualFIle, FlowVisualFIle);
            // my_flow.FlowData.PressureInteriorEdge = my_flow.FlowData.PressureInteriorEdge.setZero();
            // my_flow.FlowData.PressureEles = my_flow.FlowData.PressureEles.setZero();
            // my_flow.FlowData.VelocityNormalScalarSepEdges = my_flow.FlowData.VelocityNormalScalarSepEdges.setZero();
            // diffusion no velocity

            my_mesh.StoreInH5(MeshH5);
            my_dfn.StoreInH5(DFNH5);
            my_flow.StoreInH5(FlowH5);
            return 0;
        }
        else
        {
            my_dfn.LoadClassFromH5(DFNH5);
            my_mesh.LoadClassFromH5(MeshH5);
            my_flow.LoadClassFromH5(FlowH5);   
        };

        my_PT.PTData.IfPureDiffusion = 1;  // > 0 means pure diffusion
        my_PT.IfPeriodic = true;
        my_PT.LoadParametersFromCSV(csvPT);
        my_PT.ParticleTracking(my_dfn,
            my_mesh,
            my_flow);
        std::string PTVisualFIle = "PT_Visual";
        my_PT.Visualization(my_dfn,
            my_mesh,
            my_flow, PTVisualFIle, PTVisualFIle, FlowVisualFIle);

        return 0;
    }
    catch (cuDFNsys::ExceptionsPause &e)
    {
        std::cout << e.what() << "\n";
        
        return 0;    
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