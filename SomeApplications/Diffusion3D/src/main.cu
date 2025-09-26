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

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>

using namespace std;
namespace fs = std::filesystem;
bool IfAFileExist(const string &FilePath);
std::vector<double> readFirstLineCSV(const std::string& filename);

int main(int argc, char *argv[])
{
    if (argc != 8)
    {
        std::cout << "Usage: " << argv[0] << 
        " <csvDFN> <csvMesh> <csvFlow> <csvPT> <csvControlPlane> <IfDiscontinueAfterFirstAbsorption> <IfInletReflective>\n";
        exit(0);
    }
    std::string csvDFN = std::string(argv[1]);
    std::string csvMesh = std::string(argv[2]);
    std::string csvFlow = std::string(argv[3]);
    std::string csvPT = std::string(argv[4]);
    std::string csvControlPlanes = std::string(argv[5]);
    size_t IfDiscontinueAfterFirstAbsorption_ = atoi(argv[6]);
    size_t IfInletReflective_ = atoi(argv[7]);    

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
        bool IFDFNFILEEXIST = IfAFileExist("DFNFile.txt");

        if (!IfAFileExist(DFNH5 + ".h5") && !IFDFNFILEEXIST)
        {
            my_dfn.LoadDFNFromCSV(csvDFN);
            my_dfn.IdentifyIntersectionsClusters(true);
            if (my_dfn.PercolationCluster.size() == 0)
                return 0;
            std::cout << "It is percolative\n";
            std::string DFNVisualFIle = "DFN_Visual";
            my_dfn.Visualization(DFNVisualFIle, DFNVisualFIle, DFNVisualFIle, false, false, true, true);
            

            my_mesh.LoadParametersFromCSV(csvMesh);
            my_mesh.MeshGeneration(my_dfn);
            std::string MeshVisualFIle = "Mesh_Visual";
            my_mesh.Visualization(my_dfn, MeshVisualFIle, MeshVisualFIle, MeshVisualFIle, true, true);
            

            my_flow.LoadParametersFromCSV(csvFlow);
            my_flow.DoesNotNeedToSolve = true;
            my_flow.FlowSimulation(my_dfn, my_mesh);
            
            my_flow.Visualization(my_dfn,my_mesh, FlowVisualFIle, FlowVisualFIle, FlowVisualFIle);
            // my_flow.FlowData.PressureInteriorEdge = my_flow.FlowData.PressureInteriorEdge.setZero();
            // my_flow.FlowData.PressureEles = my_flow.FlowData.PressureEles.setZero();
            // my_flow.FlowData.VelocityNormalScalarSepEdges = my_flow.FlowData.VelocityNormalScalarSepEdges.setZero();
            // diffusion no velocity
            my_dfn.StoreInH5(DFNH5);
            my_mesh.StoreInH5(MeshH5);
            my_flow.StoreInH5(FlowH5);
            return 0;
        }
        else if (IfAFileExist(DFNH5 + ".h5") && !IFDFNFILEEXIST)
        {
            my_dfn.LoadClassFromH5(DFNH5);
            my_mesh.LoadClassFromH5(MeshH5);
            my_flow.LoadClassFromH5(FlowH5); 
            my_flow.FlowData.DoesNotNeedToSolve = true;  
        }
        else if (IFDFNFILEEXIST)
        {
            std::ifstream file("DFNFile.txt");  // open the file
            if (!file.is_open()) {
                std::cerr << "Could not open file `DFNFile.txt`!" << std::endl;
                return 1;
            }

            std::string firstLine;
            if (std::getline(file, firstLine)) {
                std::cout << "First line: " << firstLine << std::endl;
            } else {
                std::cerr << "File is empty or error reading." << std::endl;
            }

            file.close();

            my_dfn.LoadClassFromH5(firstLine + "/" + DFNH5);
            my_mesh.LoadClassFromH5(firstLine + "/" + MeshH5);
            my_flow.LoadClassFromH5(firstLine + "/" + FlowH5); 
            my_flow.FlowData.DoesNotNeedToSolve = true; 
        }
        else
        {
            std::cout << "Wrong mode of loading dfn\n";
            exit(0);
        };

        if (IfAFileExist("./ParticlePositionResult/DispersionInfo.h5") && IfDiscontinueAfterFirstAbsorption_)
        {
            cuDFNsys::HDF5API hdf5_rw;
            std::vector<uint> tmp = 
                hdf5_rw.ReadDataset<uint>("./ParticlePositionResult/DispersionInfo.h5", "N", "NumParticlesLeftFromInlet");
            if (tmp[0] > 1)
            {
                std:: cout << "Found first hit on inlet\n";
                return;
            }
            tmp = hdf5_rw.ReadDataset<uint>(
                            "./ParticlePositionResult/DispersionInfo.h5",
                        "N", "NumParticles");
            int NumParticlesTotal = tmp[0];
            std::vector<uint> FPT = hdf5_rw.ReadDataset<uint>(
                        "./ParticlePositionResult/"
                            "ParticlePosition_FPTControlPlanes.h5",
                        "N",
                        "ControlPlane_" +
                            std::to_string(my_dfn.DomainSizeX * my_dfn.DomainDimensionRatio.z * -0.5) +
                            "_m");
            int zeroCount = std::count(FPT.begin(), FPT.end(), 0);
            if (zeroCount < NumParticlesTotal)
            {
                std:: cout << "Found first hit on outlet\n";
                return;
            }
        };//
        my_PT.LoadParametersFromCSV(csvPT);

        my_PT.IfPureDiffusion = 1;  // > 0 means pure diffusion
        my_PT.IfDiscontinueAfterFirstAbsorption = IfDiscontinueAfterFirstAbsorption_; // > 0 means stop when there is first absorption
        my_PT.IfReflectionAtInlet = IfInletReflective_; // > 0 means reflective inlet

        my_PT.ControlPlanes = readFirstLineCSV(csvControlPlanes);
        // cout << "b\n";
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
        
        // int result = system("rm -rf ./*.h5");

        return 0;    
    }
    catch (cuDFNsys::ExceptionsIgnore &e)
    {
        std::cout << e.what() << "\n";
        // int result = system("rm -rf ./*.h5");
        return 0;    
    }
    catch (H5::Exception &e)
    {
        //std::cout << e.what() << "\n";
        std::cout << "H5::Exceptions: " << e.getDetailMsg() << endl;
        // int result = system("rm -rf ./*.h5");
        return 0;    
    }
    catch (H5::FileIException &e)
    {
        //std::cout << e.what() << "\n";
        std::cout << "H5::Exceptions: " << e.getDetailMsg() << endl;
        // int result = system("rm -rf ./*.h5");
        return 0;    
    }
    catch (...)
    {
        std::cout << "Unclear exception\n";
        // int result = system("rm -rf ./*.h5");
        throw;
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

std::vector<double> readFirstLineCSV(const std::string& filename1) 
{
    std::string filename = filename1 + ".csv";

    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    std::string firstLine;
    if (!std::getline(file, firstLine)) {
        throw std::runtime_error("File is empty: " + filename);
    }
    
    std::vector<double> result;
    std::istringstream lineStream(firstLine);
    std::string token;
    
    while (std::getline(lineStream, token, ',')) {
        // Trim whitespace from the token
        token.erase(0, token.find_first_not_of(" \t\n\r\f\v"));
        token.erase(token.find_last_not_of(" \t\n\r\f\v") + 1);
        
        // Stop reading if token is empty
        if (token.empty()) {
            break;
        }
        
        // Try to convert token to double
        try {
            double value = std::stod(token);
            result.push_back(value);
        }
        catch (const std::invalid_argument&) {
            // Stop reading if conversion fails (non-numeric value)
            break;
        }
        catch (const std::out_of_range&) {
            // Handle out of range values if needed, or break
            break;
        }
    }
    
    return result;
}
