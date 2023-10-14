/****************************************************************************
* cuDFNsys - simulating flow and transport in 3D fracture networks          *
* Copyright (C) 2022, Tingchang YIN, Sergio GALINDO-TORRES                  *
*                                                                           *
* This program is free software: you can redistribute it and/or modify      *
* it under the terms of the GNU Affero General Public License as            *
* published by the Free Software Foundation, either version 3 of the        *
* License, or (at your option) any later version.                           *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU Affero General Public License for more details.                       *
*                                                                           *
* You should have received a copy of the GNU Affero General Public License  *
* along with this program.  If not, see <https://www.gnu.org/licenses/>.    *
*****************************************************************************/

// ====================================================
// NAME:        A quickstart example to generate DFNs
// DESCRIPTION: Call cuDFNsys functions to do simulation.
// AUTHOR:      Tingchang YIN
// DATE:        13/10/2023
// ====================================================

#include "cuDFNsys.cuh"
#include <fstream>
#include <iostream>
#include <limits.h>
#include <unistd.h>

int main(int argc, char *argv[])
{

    int dev = 0;
    // No. 0 GPU card
    GPUErrCheck(cudaSetDevice(dev));
    // try to use the first GPU card and check its availability

    time_t t;
    time(&t);
    // t is a random seed

    // here Two sets of fractures are generated
    // Set 1: kappa = 10, mean orientation is (0, 0, 1)
    //        fracture sizes is a power law, alpha = 1.5, minR = 1, maxR = 15
    //        Number of Fractures = 100
    //        for aperture b, beta = 0.25, gamma = 5.1e-4
    //        The formula is b = gamma * R ^ beta
    // Set 2: kappa = 5, mean orientation is (1, 0, 0)
    //        fracture sizes is a lognormal distribution,
    //        log mean = 8.5, log std. dev = 5.5, minR = 1, maxR = 15
    //        Number of Fractures = 100
    //        for aperture b, beta = 0.15, gamma = 7.1e-6
    // The domain size is: X_L = 30, Y_L = 30, Z_L = 60
    // The pre-defined percolation direction is along the z axis
    // means that we check percolation state along z
    // the hydraulic head gradient along z is one
    // no flux can go out from lateral surfaces
    // now, let's define some variables to represent these parameters

    std::vector<int> NumFracture = {100, 100};
    // Number of fractures for the two sets
    std::vector<double> Kappa = {10, 5};
    // kappa values for the two sets
    std::vector<cuDFNsys::Vector3<double>> MeanOrientationOfFisherDistribution =
        {
            cuDFNsys::MakeVector3(0., 0., 1.),
            cuDFNsys::MakeVector3(1., 0., 0.)};
    // mean orientations of Fisher distribution for the two sets
    // cuDFNsys::Vector3<double> is a vector with 3 elements (i.e., x, y, z)
    double DomainSize_X = 30;
    // Domain size in the x direction
    double3 DomainDimensionRatio = make_double3(1, 1, 2);
    // Domain dimension ratio (30/30, 30/30, 60/30)
    std::vector<double> Beta = {0.25, 0.15};
    std::vector<double> Gamma = {5.1e-4, 7.1e-6};
    // Beta and Gamma for calculating the aperture of fractures
    // and conductivity is (b^3)/12
    std::vector<int> Sign_FractureSizeDistributionPattern = {0, 1};
    // clarify the distribution pattern for the two sets
    // 0 is a sign of power-laws
    // 1 denotes lognotmal distributions
    std::vector<cuDFNsys::Vector4<double>> SizeDistributionParameters =
        {
            cuDFNsys::MakeVector4(1.5, 1., 15., 0.),
            cuDFNsys::MakeVector4(8.5, 5.5, 1., 15.)};
    // parameters for the size distributions, power law and lognormal
    // cuDFNsys::Vector4<double> is a vector with 4 elements (i.e., x, y, z and w)
    int Perco_dir = 2;
    // the pre-defined percolation direction is along z axis
    // or say, flow direction is along z axis

    int NumFracturesTotal = 0;
    // the total number of fractures
    for (int i = 0; i < NumFracture.size(); ++i)
        NumFracturesTotal += NumFracture[i];

    // now, let's create empty fracture vectors with thrust vectors

    thrust::host_vector<cuDFNsys::Fracture<double>> Frac_host(NumFracturesTotal);
    // this is the vector on the host side
    thrust::device_vector<cuDFNsys::Fracture<double>> Frac_device(NumFracturesTotal);
    // this is the vector on the device side

    // cuDFNsys::Fracture<double> is a struct variable representing fractures

    cuDFNsys::Fracture<double> *Frac_device_ptr;
    Frac_device_ptr = thrust::raw_pointer_cast(Frac_device.data());
    // a pointer to device fracture vector

    double istart = cuDFNsys::CPUSecond();
    // we can count time

    // let's create fractures
    for (int i = 0; i < NumFracture.size(); ++i)
    {
        // since there are two sets of fractures
        // we generate one by one
        // generate fractures in local vectors, then copy them to the global vectors

        thrust::host_vector<cuDFNsys::Fracture<double>> Frac_host_sub(NumFracture[i]);
        thrust::device_vector<cuDFNsys::Fracture<double>> Frac_device_sub(NumFracture[i]);
        cuDFNsys::Fracture<double> *Frac_device_sub_ptr;
        Frac_device_sub_ptr = thrust::raw_pointer_cast(Frac_device_sub.data());

        cuDFNsys::Fractures<double><<<NumFracture[i] / 256 + 1, 256>>>(Frac_device_sub_ptr,                     // the pointer to device vector
                                                                       (unsigned long)t,                        // seed
                                                                       NumFracture[i],                          // number of fracture
                                                                       DomainSize_X,                            // domain size
                                                                       Sign_FractureSizeDistributionPattern[i], // distribution pattern of fracture sizes
                                                                       SizeDistributionParameters[i],           // parameters of distribution of fracture sizes
                                                                       Kappa[i],                                // kappa value of fisher distribution
                                                                       Beta[i],                                 // beta
                                                                       Gamma[i],                                // gamma
                                                                       DomainDimensionRatio,                    // ratio of domain dimensions
                                                                       MeanOrientationOfFisherDistribution[i]); // mean orientations

        cudaDeviceSynchronize();
        // wait until the device function finish
        Frac_host_sub = Frac_device_sub;

        thrust::copy(thrust::host, Frac_host_sub.begin(), Frac_host_sub.end(), Frac_host.begin() + i * NumFracture[i]);
    }

    Frac_device = Frac_host;
    // let device vector = host vector
    Frac_device_ptr = thrust::raw_pointer_cast(Frac_device.data());
    // update the pointer of device vector

    double ielaps = cuDFNsys::CPUSecond() - istart;
    cout << "Running time of fracture generation: " << ielaps << " sec\n";

    std::map<pair<size_t, size_t>, pair<cuDFNsys::Vector3<double>, cuDFNsys::Vector3<double>>> Intersection_map;
    // a map for storing intersections,
    // Intersection_map.first = a pair of fracture No.
    // Intersection_map.second = a pair of coordinates, i.e. the ends of the intersection

    istart = cuDFNsys::CPUSecond();
    cuDFNsys::IdentifyIntersection<double> identifyInters{(size_t)NumFracturesTotal, // number of fractures
                                                          Frac_device_ptr,           // pointer of device vector of fractures
                                                          false,                     // if you want to use truncated fractures? here is false,
                                                          Intersection_map};
    // this process is based on GPU
    ielaps = cuDFNsys::CPUSecond() - istart;
    cout << "Running time of identification of intersections: " << ielaps << " sec\n";

    std::vector<std::vector<size_t>> ListClusters; // will store all fracture clusters
    std::vector<size_t> Percolation_cluster;       // will store the No. of percolating cluster

    istart = cuDFNsys::CPUSecond();
    cuDFNsys::Graph<double> G{(size_t)NumFracturesTotal, Intersection_map};
    G.UseDFS(ListClusters);
    // DFS algorithm to identify clusters

    cuDFNsys::IdentifyPercolationCluster<double> IdentiClu{ListClusters,         // all clusters
                                                           Frac_host,            // host vector of fractures
                                                           Perco_dir,            // percolation direction / flow direction
                                                           Percolation_cluster}; // percolation cluster
    ielaps = cuDFNsys::CPUSecond() - istart;
    cout << "Running time of identification of clusters: " << ielaps << " sec\n";

    // now a DFN is randomly generated, let's see it
    cuDFNsys::MatlabPlotDFN<double> As{"DFN.h5",            // the .h5 file, with suffix
                                       "DFN.m",             // matlab script to visualize the DFN, with suffix
                                       Frac_host,           // host vector of fractures
                                       Intersection_map,    // intersection map
                                       ListClusters,        // clusters
                                       Percolation_cluster, // No. or say, ID, of percolating cluster
                                       false,               // if show truncated fractures?
                                       true,                // if show intersections?
                                       true,                // if show clusters?
                                       true,                // if show orientations data?
                                       DomainSize_X,        // domain size
                                       Perco_dir,           // flow direction
                                       true,                // true means I also want to see DFN with python script, a .py file will be generated
                                       "DFN",               // the name of python script, without suffix
                                       DomainDimensionRatio};

    // then I want to identify intersections of TRUNCATED fractures!!!
    // as well as the percolation cluster
    // this is important for flow simulation!!!

    Intersection_map.clear();
    ListClusters.clear();
    Percolation_cluster.clear();
    // for the sake of safety, clear them first

    /// now let's consider truncated fractures!!!
    istart = cuDFNsys::CPUSecond();
    cuDFNsys::IdentifyIntersection<double> identifyInters2___{(size_t)NumFracturesTotal,
                                                              Frac_device_ptr,
                                                              true,
                                                              Intersection_map};
    ielaps = cuDFNsys::CPUSecond() - istart;
    cout << "Running time of identification of intersections with truncated fractures: " << ielaps << " sec\n";

    istart = cuDFNsys::CPUSecond();
    cuDFNsys::Graph<double> G2__{(size_t)NumFracturesTotal, Intersection_map};
    G2__.UseDFS(ListClusters);
    // DFS algorithm to identify clusters

    cuDFNsys::IdentifyPercolationCluster<double> IdentiClu2____{ListClusters,
                                                                Frac_host,
                                                                Perco_dir,
                                                                Percolation_cluster};
    ielaps = cuDFNsys::CPUSecond() - istart;
    cout << "Running time of identification of clusters with truncated fractures: " << ielaps << " sec\n";

    cuDFNsys::MatlabPlotDFN<double> As2__{"DFN_truncated.h5",  // the .h5 file, with suffix
                                          "DFN_truncated.m",   // matlab script to visualize the DFN, with suffix
                                          Frac_host,           // host vector of fractures
                                          Intersection_map,    // intersection map
                                          ListClusters,        // clusters
                                          Percolation_cluster, // No. or say, ID, of percolating cluster
                                          true,                // if show truncated fractures?
                                          true,                // if show intersections?
                                          true,                // if show clusters?
                                          true,                // if show orientations data?
                                          DomainSize_X,        // domain size
                                          Perco_dir,           // flow direction
                                          true,                // true means I also want to see DFN with python script, a .py file will be generated
                                          "DFN_truncated",     // the name of python script, without suffix
                                          DomainDimensionRatio};

    if (Percolation_cluster.size() > 0) // please be sure that there is at least one spanning cluster
    {
        // then output the fracture information for next simulation (mesh, flow, transport)

        cuDFNsys::OutputObjectData<double> lk_out;
        lk_out.OutputFractures("Fractures.h5", Frac_host, DomainSize_X, DomainDimensionRatio);
        cout << "***** The DFN is percolative, which has been output *****\n";
    };

    return 0;
};