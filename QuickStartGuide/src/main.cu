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
// NAME:        A quickstart example
// DESCRIPTION: Call cuDFNsys functions to do simulation.
// AUTHOR:      Tingchang YIN
// DATE:        03/11/2022
// ====================================================

#include "cuDFNsys.cuh"
#include <unistd.h>

#ifdef USE_DOUBLES
typedef double _DataType_;
#else
typedef float _DataType_;
#endif
// _DataType_ here is actually double

int main(int argc, char *argv[])
{

    try
    {
        int dev = 0;                     // No. 0 GPU card
        GPUErrCheck(cudaSetDevice(dev)); // try to use the first GPU card and check its availability

        int DSIZE = 150;
        // the number of fractures is 150

        _DataType_ L = 30;
        // the domain size is 30 m, the center of the domain is (0, 0, 0)

        _DataType_ minGrid = 2;
        _DataType_ maxGrid = 3;
        // we set the minimum grid size of mesh is 2m, maxmum is 3 m.
        // they are just input parameter, hence the actual grid size may differ from them more or less

        _DataType_ kappa_para = 0;
        // here we just set fisher constant = 0
        // i.e. uniform orientation

        _DataType_ beta = 0.25;
        // the formula is: b = gamma * r ^ beta, b is aperture
        // r is fracture size
        // gamma here is 5e-4, a default value in the function cuDFNsys::Fractures

        int ModeSizeDistri = 0;
        // ModeSizeDistri here means that the distribution of fracture sizes is a power law
        // fracture size here is the radius of fractures

        cuDFNsys::Vector4<_DataType_> ParaSizeDistri =
            cuDFNsys::MakeVector4(1.5,
                                  1.,
                                  15.,
                                  0.);
        // ParaSizeDistri defines the parameter of distribution of fracture sizes
        // 1.5 = the exponent of a power law
        // 1. = the minimum fracture size
        // 15. = the maximum

        int perco_dir = 2;
        // the pre-defined percolation direction is along z axis
        // or say, flow direction is along z axis

        thrust::host_vector<cuDFNsys::Fracture<_DataType_>> Frac_verts_host(DSIZE);
        thrust::device_vector<cuDFNsys::Fracture<_DataType_>> Frac_verts_device(DSIZE);
        cuDFNsys::Fracture<_DataType_> *Frac_verts_device_ptr;
        Frac_verts_device_ptr = thrust::raw_pointer_cast(Frac_verts_device.data());
        // create host and device vectors of cuDFNsys::Fracture
        // here we use double precision, namely, _DataType_
        // a device pointer, 'Frac_verts_device_ptr', pointing to the elements in device vector

        cuDFNsys::Warmup<<<DSIZE / 256 + 1, 256>>>();
        cudaDeviceSynchronize();
        // let us warmup the GPU firstly
        // cudaDeviceSynchronize() is very important, means that,
        // after the above kernel function is finished, then the following scripts are implemented.

        time_t t;
        time(&t);
        // t is for random seed

        double istart = cuDFNsys::CPUSecond(); // we can count time

        cuDFNsys::Fractures<_DataType_><<<DSIZE / 256 + 1, 256>>>(Frac_verts_device_ptr, // the pointer to device vector
                                                                  (unsigned long)t,      // seed
                                                                  DSIZE,                 // number of fracture
                                                                  L,                     // domain size
                                                                  ModeSizeDistri,        // distribution pattern of fracture sizes
                                                                  ParaSizeDistri,        // parameters of distribution of fracture sizes
                                                                  kappa_para,            // kappa value of fisher distribution
                                                                  beta);                 // beta

        cudaDeviceSynchronize();             // now we finished generating fractures
        Frac_verts_host = Frac_verts_device; // copy data from device to host

        double ielaps = cuDFNsys::CPUSecond() - istart;
        cout << "Running time of fracture generation: " << ielaps << " sec\n";

        std::map<pair<size_t, size_t>, pair<cuDFNsys::Vector3<_DataType_>, cuDFNsys::Vector3<_DataType_>>> Intersection_map;
        // a map for storing intersections,
        // Intersection_map.first = a pair of fracture No.
        // Intersection_map.second = a pair of coordinates, i.e. the ends of the intersection

        istart = cuDFNsys::CPUSecond();
        cuDFNsys::IdentifyIntersection<_DataType_> identifyInters{Frac_verts_host.size(), // number of fractures
                                                                  Frac_verts_device_ptr,  // pointer of device vector of fractures
                                                                  false,                  // if you want to use truncated fractures? here is false,
                                                                  Intersection_map};
        // this process is based on GPU
        ielaps = cuDFNsys::CPUSecond() - istart;
        cout << "Running time of identification of intersections: " << ielaps << " sec\n";

        std::vector<std::vector<size_t>> ListClusters; // will store all fracture clusters
        std::vector<size_t> Percolation_cluster;       // will store the No. of percolating cluster

        istart = cuDFNsys::CPUSecond();
        cuDFNsys::Graph<_DataType_> G{(size_t)DSIZE, Intersection_map};
        G.UseDFS(ListClusters);
        // DFS algorithm to identify clusters

        cuDFNsys::IdentifyPercolationCluster<_DataType_> IdentiClu{ListClusters,         // all clusters
                                                                   Frac_verts_host,      // host vector of fractures
                                                                   perco_dir,            // percolation direction / flow direction
                                                                   Percolation_cluster}; // percolation cluster
        ielaps = cuDFNsys::CPUSecond() - istart;
        cout << "Running time of identification of clusters: " << ielaps << " sec\n";

        // now a DFN is randomly generated, let's see it
        cuDFNsys::MatlabPlotDFN<_DataType_> As{"DFN_I.h5",          // the .h5 file, with suffix
                                               "DFN_I.m",           // matlab script to visualize the DFN, with suffix
                                               Frac_verts_host,     // host vector of fractures
                                               Intersection_map,    // intersection map
                                               ListClusters,        // clusters
                                               Percolation_cluster, // No. or say, ID, of percolating cluster
                                               false,               // if show truncated fractures?
                                               true,                // if show intersections?
                                               true,                // if show clusters?
                                               true,                // if show orientations data?
                                               L,                   // domain size
                                               perco_dir,           // flow direction
                                               true,                // true means I also want to see DFN with python script, a .py file will be generated
                                               "DFN_I"};            // the name of python script, without suffix

        // then I want to identify intersections of TRUNCATED fractures!!!
        // as well as the percolation cluster
        // this is important for flow simulation!!!

        Intersection_map.clear();
        ListClusters.clear();
        Percolation_cluster.clear();
        // for the sake of safety, clear them first

        /// now let's consider truncated fractures!!!
        istart = cuDFNsys::CPUSecond();
        cuDFNsys::IdentifyIntersection<_DataType_> identifyInters2___{Frac_verts_host.size(),
                                                                      Frac_verts_device_ptr,
                                                                      true,
                                                                      Intersection_map};
        ielaps = cuDFNsys::CPUSecond() - istart;
        cout << "Running time of identification of intersections with truncated fractures: " << ielaps << " sec\n";

        istart = cuDFNsys::CPUSecond();
        cuDFNsys::Graph<_DataType_> G2__{(size_t)DSIZE, Intersection_map};
        G2__.UseDFS(ListClusters);
        // DFS algorithm to identify clusters

        cuDFNsys::IdentifyPercolationCluster<_DataType_> IdentiClu2____{ListClusters,
                                                                        Frac_verts_host,
                                                                        perco_dir,
                                                                        Percolation_cluster};
        ielaps = cuDFNsys::CPUSecond() - istart;
        cout << "Running time of identification of clusters with truncated fractures: " << ielaps << " sec\n";

        cuDFNsys::MatlabPlotDFN<_DataType_> As2__{"DFN_II.h5",         // the .h5 file, with suffix
                                                  "DFN_II.m",          // matlab script to visualize the DFN, with suffix
                                                  Frac_verts_host,     // host vector of fractures
                                                  Intersection_map,    // intersection map
                                                  ListClusters,        // clusters
                                                  Percolation_cluster, // No. or say, ID, of percolating cluster
                                                  true,                // if show truncated fractures?
                                                  true,                // if show intersections?
                                                  true,                // if show clusters?
                                                  true,                // if show orientations data?
                                                  L,                   // domain size
                                                  perco_dir,           // flow direction
                                                  true,                // true means I also want to see DFN with python script, a .py file will be generated
                                                  "DFN_II"};           // the name of python script, without suffix

        Frac_verts_device.clear();
        Frac_verts_device.shrink_to_fit();
        // device vector of fractures now is not neccessary

        if (Percolation_cluster.size() > 0) // please be sure that there is at least one spanning cluster
        {
            std::vector<size_t> Fracs_percol; // will store ID / No of fractures in percolating cluster

            istart = cuDFNsys::CPUSecond();
            cuDFNsys::GetAllPercolatingFractures GetPer{Percolation_cluster,
                                                        ListClusters,
                                                        Fracs_percol};
            // this function is simple, just collecting fractures in the percolating clusters

            std::vector<pair<int, int>> IntersectionPair_percol;
            // will store the intersection pair of fractures in percolation clusters

            cuDFNsys::RemoveDeadEndFrac<_DataType_> RDEF{Fracs_percol,            // fractures' ID in percolating cluster
                                                         IntersectionPair_percol, // intersection pair
                                                         (size_t)perco_dir,       // flow direction
                                                         Frac_verts_host,         // host vector of fractures
                                                         Intersection_map};       // map of intersection
            // the above function removes dead end fractures

            ielaps = cuDFNsys::CPUSecond() - istart;
            cout << "Running time of removing dead-end fractures: " << ielaps << " sec\n";

            istart = cuDFNsys::CPUSecond();
            cuDFNsys::Mesh<_DataType_> mesh{Frac_verts_host,         // host vector of fractures, after removing fractures of dead end
                                            IntersectionPair_percol, // intersection pair
                                            &Fracs_percol,           // fractures' ID in percolating cluster
                                            minGrid,                 // minimum grid size
                                            maxGrid,                 // maximum grid size
                                            perco_dir,               // flow direction
                                            L};                      // domain size
            // mesh finished
            mesh.MatlabPlot("DFN_mesh.h5",   // h5 file
                            "DFN_mesh.m",    // name of matlab script
                            Frac_verts_host, // fracture vector on the host side
                            L,               // domain size
                            true,            // if check 2D coordinates, because 3D fractures can be mapped to 2D plane
                            true,            // if check the edge attributes? Neumann, Dirichlet?
                            true,            // if I want to see mesh with Python?
                            "DFN_mesh");     // name of python script without suffix
            ielaps = cuDFNsys::CPUSecond() - istart;
            cout << "Running time of mesh: " << ielaps << " sec\n";

            istart = cuDFNsys::CPUSecond();
            cuDFNsys::MHFEM<_DataType_> fem{mesh,            // mesh object
                                            Frac_verts_host, // fractures
                                            100,             // hydraulic head of inlet = 100 m
                                            20,              // hydraulic head of outlet = 20 m
                                            perco_dir,       // flow direction
                                            L};              // domain size

            fem.MatlabPlot("MHFEM.h5",      // h5 file
                           "MHFEM.m",       // matlab script to see mhfem result
                           Frac_verts_host, // fractures
                           mesh,            // mesh object
                           L,               // domain size
                           true,            // if use python to do visualization
                           "MHFEM");        // name of python script, without suffix

            ielaps = cuDFNsys::CPUSecond() - istart;
            cout << "Running time of mhfem: " << ielaps << " sec\n";

            cuDFNsys::OutputObjectData<_DataType_> lk;
            lk.OutputFractures("FracturesForParticle.h5", Frac_verts_host, L);
            // the above two command just in order to output fractures information to transform 2D particle data to 3D

            istart = cuDFNsys::CPUSecond();
            cuDFNsys::ParticleTransport<_DataType_> p{(unsigned long)t,          // random seed
                                                      atoi(argv[1]),             // number of particle
                                                      atoi(argv[2]),             // number of time steps
                                                      (_DataType_)atof(argv[3]), // delta T
                                                      (_DataType_)atof(argv[4]), // molecular diffusion
                                                      Frac_verts_host,           // fractures
                                                      mesh,                      // mesh object
                                                      fem,                       // mhfem object
                                                      (uint)perco_dir,           // flow direction
                                                      -0.5 * L,                  // the target plane, z = -0.5 * L
                                                      "Particle_tracking",       // use particle tracking algorithm
                                                      "Flux-weighted"};          // the injection mode is flux-weighted

            p.MatlabPlot("MHFEM.h5",            // h5 file of mhfem
                         "ParticlesMovement.m", // matlab script
                         mesh,                  // mesh result
                         fem,                   // mhfem object
                         L,                     // domain size
                         true,
                         "ParticlesMovement");

            // note that right now data in the output file () is not 3D
            // you have to transform 2D data to 3D
            // which can be done by run the executable file 'Transform2DH5ParticleDataTo3D'
            // then, the ParticlesMovement.m can be implemented in matlab
            // or ParticlesMovement.py can be run
            // The 'compileCode.sh' will compile the 'Transform2DH5ParticleDataTo3D' code and run it after data (particle positions) output

            ielaps = cuDFNsys::CPUSecond() - istart;
            cout << "Running time of particle tracking: " << ielaps << " sec\n";
        }
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
    };
    return 0;
};