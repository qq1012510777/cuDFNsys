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

template<typename T>
bool isInVector(const std::vector<T>& vec, const T& element) {
    return std::find(vec.begin(), vec.end(), element) != vec.end();
}

int main(int argc, char *argv[])
{
    time_t t;
    time(&t);
    if (argc != 5)
    {
        std::cout << "usgae: " << argv[0] << " <domain size> <csv name without suffix> <min grid size> <max grid size>\n";
        exit(0);
    }
    double sizeDomain = atof(argv[1]);
    cuDFNsys::DFN<double> my_dfn;
    my_dfn.LoadDFNFromCSV(std::string(argv[2]));
    my_dfn.IdentifyIntersectionsClusters(false);
    my_dfn.StoreInH5("DFN_showCritical_"+std::to_string(my_dfn.RandomSeed));
    my_dfn.Visualization("DFNVisual", "DFNVisual", "DFNVisual", false, false, true, true);

    cuDFNsys::MeshDFN<double> mesh;
    mesh.MinElementSize = atof(argv[3]);
    mesh.MaxElementSize = atof(argv[4]);

    mesh.MeshGeneration(my_dfn);
    mesh.Visualization(my_dfn, "Mesh_show", "Mesh_show", "Mesh_show", false, false);

    cuDFNsys::FlowDFN<double> flow;
    flow.InletHead = 1;
    flow.OutletHead = 0;
    flow.FlowSimulation(my_dfn,
        mesh);
    flow.Visualization(my_dfn, mesh, "Flow_show", "Flow_show", "Flow_show");
    std::cout << "flow.FlowData.Permeability=" << flow.FlowData.Permeability << "\n";
    return 0;
}
