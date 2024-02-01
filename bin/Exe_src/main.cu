#include "cuDFNsys.cuh"
#include <filesystem>
string ColoringString(const string &s);
string PrintAllCsvAndSelectone();
int main()
{
    try
    {
        time_t t;
        time(&t);

        cout << "\n\n"
             << ColoringString("************** cuDFNsys starts **************")
             << endl;
        uint options = 0;

        std::ifstream File_d("./Class_DFN.h5");
        bool DFSW = File_d.good();

        cuDFNsys::DFN<double> my_dfn;
        my_dfn.RandomSeed = (unsigned long)t;

        if (DFSW)
        {
            cout << ColoringString(
                        "A DFN file is existing. Do you want to load this DFN?")
                 << endl;
            cout << ColoringString("Input 0 for No, 1 for yes:") << endl;
            cin >> options;
            if (options == 0)
                goto NewDFN;
            else
            {
                cout << ColoringString("Loading an existing DFN") << endl;
                my_dfn.LoadClassFromH5("Class_DFN");
                my_dfn.IdentifyIntersectionsClusters(true);
            }
        }
        else
        {
        NewDFN:;
            cout << ColoringString("Creating a new DFN. Please provide the "
                                   "name of the .csv\n\n")
                 << endl;
            string nameCSV = PrintAllCsvAndSelectone();
            //cin >> nameCSV;

            my_dfn.LoadDFNFromCSV(nameCSV);
            my_dfn.IdentifyIntersectionsClusters(true);
        }
        File_d.close();
        //my_dfn.StoreInH5("Class_DFN");
        my_dfn.Visualization("DFN_VISUAL", "DFN_VISUAL", "DFN_VISUAL", true,
                             true, true, true);

        //---------------------mesh------------
        cuDFNsys::MeshDFN<double> meshGen;

        std::ifstream File_d2("./Class_MESH.h5");
        DFSW = File_d2.good();

        if (DFSW)
        {
            cout << ColoringString("A mesh file is existing! Are you sure that "
                                   "it matchs the "
                                   "DFN?")
                 << endl;
            cout << ColoringString("Input 0 for No, 1 for yes:") << endl;
            cin >> options;
            if (options == 0)
                goto NewMesh;
            else
            {
                cout << ColoringString("Loading an existing mesh") << endl;
                meshGen.LoadClassFromH5("Class_MESH");
            }
        }
        else
        {
        NewMesh:;
            cout << ColoringString("Do you want to generate mesh? Input 0 for "
                                   "No, 1 for yes:")
                 << endl;
            cin >> options;
            if (options == 0)
                return 0;

            cout << ColoringString("Creating a new mesh") << endl;
            // cout << ColoringString("Input the minimum grid size:") << endl;
            // cin >> meshGen.MinElementSize;
            // cout << ColoringString("Input the maximum grid size:") << endl;
            // cin >> meshGen.MaxElementSize;
            meshGen.LoadParametersFromCSV("Mesh_parameters");
            meshGen.MeshGeneration(my_dfn);
        }

        File_d2.close();

        meshGen.Visualization(my_dfn, "DFN_MESH_VISUAL", "DFN_MESH_VISUAL",
                              "DFN_MESH_VISUAL", true, true);

        my_dfn.StoreInH5("Class_DFN");
        meshGen.StoreInH5("Class_MESH");

        //-------------------- flow
        cuDFNsys::FlowDFN<double> flowDFN;

        std::ifstream File_d3("./Class_FLOW.h5");
        DFSW = File_d3.good();

        if (DFSW)
        {
            cout << ColoringString("A flow file is existing! Are you sure that "
                                   "it matchs the "
                                   "DFN?")
                 << endl;
            cout << ColoringString("Input 0 for No, 1 for yes:") << endl;
            cin >> options;
            if (options == 0)
                goto NewFlow;
            else
            {
                cout << ColoringString("Loading an existing flow file") << endl;
                flowDFN.LoadClassFromH5("Class_FLOW");
            }
        }
        else
        {
        NewFlow:;
            cout << ColoringString(
                        "Do you want to solve flow? Input 0 for No, 1 for yes:")
                 << endl;
            cin >> options;
            if (options == 0)
                return 0;

            cout << ColoringString("solving flow") << endl;
            flowDFN.LoadParametersFromCSV("Flow_parameters");
            flowDFN.FlowSimulation(my_dfn, meshGen);
        }

        flowDFN.Visualization(my_dfn, meshGen, "DFN_FLOW_VISUAL",
                              "DFN_FLOW_VISUAL", "DFN_FLOW_VISUAL");
        flowDFN.StoreInH5("Class_FLOW");
        File_d3.close();

        cout << ColoringString(
                    "\n\nmean fracture velocity is " +
                    cuDFNsys::ToStringWithWidth(flowDFN.MeanVelocity, 8))
             << endl;
        cout << ColoringString(
                    "max fracture velocity is " +
                    cuDFNsys::ToStringWithWidth(flowDFN.MaxVelocity, 8))
             << endl;
        cout << ColoringString(
                    "mean grid size is " +
                    cuDFNsys::ToStringWithWidth(meshGen.MeanGridSize, 8) +
                    " m^2")
             << endl;

        //------------------PT --------------

        cout << ColoringString("Do you want to run particle tracking?") << endl;
        cout << ColoringString("0 for No, 1 for yes:") << endl;
        cin >> options;
        if (options == 0)
            exit(0);

        cuDFNsys::PTDFN<double> particleTracking;
        particleTracking.LoadParametersFromCSV("PT_parameters");
        particleTracking.ParticleTracking(my_dfn, meshGen, flowDFN);
        particleTracking.Visualization(my_dfn, meshGen, flowDFN,
                                       "DFN_DISPERSION", "DFN_DISPERSION",
                                       "DFN_FLOW_VISUAL");

        cout << ColoringString("Do you want to transform the 2D particles' "
                               "spatial information to 3D. This should be done "
                               "untill you have run enough time of PT")
             << endl;
        cout << ColoringString("0 for No, 1 for yes:") << endl;
        cin >> options;
        if (options == 1)
            if (particleTracking.OutputAllPTInformationOrFPTCurve)
            {
                system("./Transform2DH5ParticleDataTo3D 0 "
                       "DFN_MESH_VISUAL.h5");
            }
            else
            {
                system("./Transform2DH5ParticleDataTo3D 1 "
                       "DFN_MESH_VISUAL.h5");
            };
    }
    catch (cuDFNsys::ExceptionsIgnore &e)
    {
        cout << e.what() << endl;
    }
    catch (cuDFNsys::ExceptionsPause &e)
    {
        cout << e.what() << endl;
    }
    catch (...)
    {
        throw;
    }

    return 0;
};

string ColoringString(const string &s) { return "\033[1;32m" + s + "\033[0m"; };
string PrintAllCsvAndSelectone()
{
    // Get the current working directory
    std::filesystem::path currentPath = std::filesystem::current_path();

    // Create a vector to store the names of CSV files
    std::vector<std::string> csvFiles;

    for (const auto &entry : std::filesystem::directory_iterator(currentPath))
        if (entry.path().extension() == ".csv")
            csvFiles.push_back(entry.path().filename().string());

    std::cout << ColoringString(".csv Files in the current directory:\n");
    for (size_t i = 0; i < csvFiles.size(); ++i)
        std::cout << "\tcsv number " << ColoringString(std::to_string(i))
                  << ": " << csvFiles[i] << "\n";

    std::cout << ColoringString(
        "Enter the number corresponding to the CSV file (must be for "
        "DFN generation) you "
        "want to select:\n");
    int userInput;
    std::cin >> userInput;

    if (userInput >= 0 && static_cast<size_t>(userInput) <= csvFiles.size())
    {
        // Get the selected CSV file
        string nameCSV = csvFiles[userInput];

        std::cout << "\tYou selected: " << nameCSV << "\n";
        size_t dotPosition = nameCSV.find_last_of('.');
        nameCSV = nameCSV.substr(0, dotPosition);

        return nameCSV;
    }
    else
    {
        std::cout << "Invalid selection.\n";
        exit(0);
    }
    return "N";
};