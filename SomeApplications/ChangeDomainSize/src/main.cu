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

int main(int argc, char *argv[])
{
    time_t t;
    time(&t);
   
    cuDFNsys::DFN<double> my_dfn;
   
    my_dfn.RandomSeed = (unsigned long)t;
    my_dfn.DomainSizeX = 60;
    my_dfn.DomainDimensionRatio = make_double3(1, 1, 1);
   
    my_dfn.NumFractures = {600};
    my_dfn.Kappa = {0};
    my_dfn.MeanOrientationOfFisherDistribution = {make_double3(0, 0, 1)};
    my_dfn.Beta = {0.2};
    my_dfn.Gamma = {1e-8};
    my_dfn.ModeOfSizeDistribution = {0};
    my_dfn.SizeDistributionParameters = {make_double4(1.5, 1, 15, 0)};
   
    my_dfn.PercoDir = 2;
    
    my_dfn.FractureGeneration();
    my_dfn.IdentifyIntersectionsClusters(false);
    my_dfn.Visualization("DFNvisual", "DFNvisual", "DFNvisual", true, true, true, true);
    my_dfn.StoreInH5("ClassDFN");

    my_dfn.ChangeDomainSize(30, make_double3(1, 1, 1));
    my_dfn.IdentifyIntersectionsClusters(false);
    my_dfn.Visualization("DFNvisualII", "DFNvisualII", "DFNvisualII", false, true, true, true);
    my_dfn.StoreInH5("ClassDFNII");
};