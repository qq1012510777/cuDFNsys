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
// NAME:        DispersionAtOneDensityValue.cu
// DESCRIPTION: Dispersion in a DFN with INCREASING percolation parameter values
// AUTHOR:      Tingchang YIN
// DATE:        02/05/2023
// ====================================================

#include "cuDFNsys.cuh"
#include <H5Exception.h>
#include <fstream>
#include <iostream>
#include <limits.h>
#include <unistd.h>
#ifdef USE_DOUBLES
typedef double _DataType_;
#else
typedef float _DataType_;
#endif

int GetRowsOfDataset(string Filename, string Datasetname);

int main(int argc, char *argv[])
{
    char cwd[PATH_MAX];
    if (getcwd(cwd, sizeof(cwd)) != NULL)
        printf("Current working dir: %s\n", cwd);
    else
        throw cuDFNsys::ExceptionsPause("getcwd() error");
    string curPath = cwd;
    //----------------------

    int dev = 0;
    GPUErrCheck(cudaSetDevice(dev));
    cuDFNsys::Warmup<<<256 / 256 + 1, 256 /*  1, 2*/>>>();
    cudaDeviceSynchronize();

    //------------Density increases
    uint InitLoop = atoi(argv[1]);
    uint FinalLoop = atoi(argv[2]);
    uint InitDensity = atoi(argv[3]);
    uint DensityIncreament = atoi(argv[4]);
    uint MaxTranLoopTimes = atoi(argv[5]);

    uint LoopTimes = FinalLoop - InitLoop + 1;

    //------------ other inputs
    _DataType_ L = atof(argv[6]);
    _DataType_ kappa_ = atof(argv[7]),
               beta_ = atof(argv[8]),
               gamma_ = atof(argv[9]);
    int size_frac_mode = atoi(argv[10]); // mode of fracture size distributions
    cuDFNsys::Vector4<_DataType_> ParaSizeDistri =
        cuDFNsys::MakeVector4((_DataType_)atof(argv[11]),
                              (_DataType_)atof(argv[12]),
                              (_DataType_)atof(argv[13]),
                              (_DataType_)atof(argv[14]));
    double3 DomainDimensionRatio = make_double3(1, 1, atof(argv[15]));
    int IfRemoveDeadEnd = atoi(argv[16]);
    _DataType_ minGridSize = atof(argv[17]);
    _DataType_ maxGridSize = atof(argv[18]);

    int NumTimeSteps_Dispersion = atoi(argv[19]);
    int NumParticlesRandomWalk = atoi(argv[20]);
    _DataType_ DeltaT = 0;
    _DataType_ Factor_mean_time_in_grid = atof(argv[21]);
    // the mean time (a characteristic grid length over the mean velocity (m/s)) for a random walk to cross a characteristic grid length
    // but this mean time was reduced, i.e., dividing by a factor (> 1)
    // then the mean time is DeltaT
    _DataType_ DiffusionLocal = 0;
    _DataType_ LengthScale_Over_Pe = 0;
    _DataType_ LengthScale = atof(argv[22]);
    _DataType_ Pe = atof(argv[23]);
    _DataType_ ControlPlaneSpacing = atof(argv[24]);
    bool IfoutputMsd = atoi(argv[25]) == 0 ? false : true;
    bool IfoutputParticleInfoAllsteps = atoi(argv[26]) == 0 ? false : true;
    int ThresholdToStop = atoi(argv[27]);
    int ThresholdForMaximumLeftParticles = atoi(argv[28]);

    string recordMode = IfoutputParticleInfoAllsteps == false ? "FPTCurve" : "OutputAll";
    _DataType_ P_in = L, P_out = 0;

    cout << "L: " << L << endl;
    cout << "Kappa: " << kappa_ << endl;
    cout << "Beta: " << beta_ << endl;
    cout << "Gamma: " << gamma_ << endl;
    cout << "Mode of fracture size distributions: " << size_frac_mode << endl;
    cout << "Parameters of the size distribution: " << ParaSizeDistri.x << ", " << ParaSizeDistri.y << ", " << ParaSizeDistri.z << ", " << ParaSizeDistri.w << endl;
    cout << "Domain's dimension ratio: " << DomainDimensionRatio.x << ", " << DomainDimensionRatio.y << ", " << DomainDimensionRatio.z << endl;
    cout << "If remove the dead ends: " << (IfRemoveDeadEnd == 0 ? "false" : "true") << endl;
    cout << "Min grid size: " << minGridSize << endl;
    cout << "Max grid size: " << maxGridSize << endl;
    cout << "Hydraulic head at the inlet and outlet: " << P_in << ", " << P_out << endl;
    cout << "Number of time steps for random walks: " << NumTimeSteps_Dispersion << endl;
    cout << "Number of particles: " << NumParticlesRandomWalk << endl;
    cout << "Factor_mean_time_in_grid: " << Factor_mean_time_in_grid << endl;
    cout << "LengthScale: " << LengthScale << endl;
    cout << "Pe: " << Pe << endl;
    cout << "The spacing of control planes: " << ControlPlaneSpacing << endl;
    cout << "IfoutputMsd: " << (IfoutputMsd == true ? "true" : "false") << endl;
    cout << "IfoutputParticleInfoAllsteps: " << (IfoutputParticleInfoAllsteps == false ? "FPTCurve" : "OutputAll") << endl;

    int perco_dir = 2;
    LengthScale_Over_Pe = LengthScale / Pe;

    //-----------------------
    cuDFNsys::HDF5API h5g;

    for (uint i = InitLoop; i <= LoopTimes + InitLoop; ++i)
    {
        string path2 = "DFN_" + cuDFNsys::ToStringWithWidth(i, 3);
        string command1 = "mkdir -p " + path2;
        system(command1.c_str());

        string command2 = curPath + "/" + path2;
        chdir(command2.c_str());

        uint DSIZE = InitDensity + (i - 1) * DensityIncreament;

        for (uint j = 0; j < MaxTranLoopTimes; j++)
        {
            try
            {
                system("echo \" \" > ../log.txt");
                cout << "DSIZE: " << DSIZE << ", j = " << j << endl;
                //-----------if a mesh h5 exists ----
                // then the DFN is percolative
                bool If_percolative = false;
                std::ifstream fileer("mesh.h5");
                If_percolative = fileer.good();

                //----------if a mesh exists
                //---------- the dispersion simultion may have been finished
                //------ let's check it
                if (If_percolative)
                {
                    std::ifstream FileYe("ParticlePositionResult/ParticlePositionLastStep.h5");
                    bool DFSW = FileYe.good();

                    if (DFSW) // that means the simulation has been conducted
                    {

                        std::vector<uint> NUMsTEPS = h5g.ReadDataset<uint>("ParticlePositionResult/DispersionInfo.h5",
                                                                           "N", "NumOfSteps");

                        int rows_ = GetRowsOfDataset("ParticlePositionResult/ParticlePositionLastStep.h5",
                                                     "Step_" + cuDFNsys::ToStringWithWidth<int>(NUMsTEPS[0], 10));
                        if (rows_ <= ThresholdToStop)
                            break; // get enough particles arrived

                        std::vector<uint> NumPaLeft = h5g.ReadDataset<uint>("ParticlePositionResult/DispersionInfo.h5",
                                                                            "N", "NumParticlesLeftFromInlet");
                        if (NumPaLeft[0] > ThresholdForMaximumLeftParticles)
                            throw cuDFNsys::ExceptionsIgnore("Too many particles left from the inlet!\n");
                    }

                    //----------amend MSD data
                    system("python ../scripts/DelDataSet.py");
                }

                //------------need more simulations
                time_t t;
                time(&t);

                srand((unsigned int)time(0));

                Eigen::MatrixXd Ter = Eigen::MatrixXd::Random(1, 1);

                thrust::host_vector<cuDFNsys::Fracture<_DataType_>> Frac_verts_host(DSIZE);
                thrust::device_vector<cuDFNsys::Fracture<_DataType_>> Frac_verts_device(DSIZE);
                cuDFNsys::Fracture<_DataType_> *Frac_verts_device_ptr;
                Frac_verts_device_ptr = thrust::raw_pointer_cast(Frac_verts_device.data());

                cuDFNsys::InputObjectData<_DataType_> lk;

                if (!If_percolative)
                {
                    cuDFNsys::Fractures<_DataType_><<<DSIZE / 256 + 1, 256>>>(Frac_verts_device_ptr,
                                                                              (unsigned long)t + (unsigned long)ceil(abs(Ter(0, 0)) * ((unsigned long)t * 1.0)),
                                                                              DSIZE, L,
                                                                              size_frac_mode,
                                                                              ParaSizeDistri,
                                                                              kappa_, // kappa
                                                                              beta_,  // beta
                                                                              gamma_, // gamma
                                                                              DomainDimensionRatio);
                    cudaDeviceSynchronize();
                    Frac_verts_host = Frac_verts_device;
                }
                else
                {
                    Frac_verts_host.resize(0);
                    lk.InputFractures("Fractures.h5", Frac_verts_host, L, DomainDimensionRatio);
                    Frac_verts_host.shrink_to_fit();
                    Frac_verts_device = Frac_verts_host;
                }

                cout << "identifying intersections with complete fractures" << endl;
                std::map<pair<size_t, size_t>, pair<cuDFNsys::Vector3<_DataType_>, cuDFNsys::Vector3<_DataType_>>> Intersection_map;
                cuDFNsys::IdentifyIntersection<_DataType_> identifyInters{Frac_verts_host.size(),
                                                                          Frac_verts_device_ptr,
                                                                          false,
                                                                          Intersection_map};
                cout << "identifying cluster with complete fractures" << endl;
                std::vector<std::vector<size_t>> ListClusters;
                std::vector<size_t> Percolation_cluster;
                cuDFNsys::Graph<_DataType_> G{(size_t)DSIZE, Intersection_map};
                G.UseDFS(ListClusters);
                cuDFNsys::IdentifyPercolationCluster<_DataType_> IdentiClu{ListClusters,
                                                                           Frac_verts_host, perco_dir,
                                                                           Percolation_cluster};
                cout << "DFN I finished" << endl;
                cuDFNsys::MatlabPlotDFN<_DataType_> As{"DFN_I.h5", "DFN_I.m",
                                                       Frac_verts_host, Intersection_map, ListClusters,
                                                       Percolation_cluster, false, true, true, true,
                                                       L, perco_dir, true, "DFN_I", DomainDimensionRatio};

                std::vector<double> Data1(14);
                cuDFNsys::GetStatistics<double>(Frac_verts_host,
                                                Intersection_map,
                                                ListClusters,
                                                Percolation_cluster, L, Data1[0], Data1[1], Data1[2], Data1[3],
                                                Data1[4], Data1[5], Data1[6], Data1[7], Data1[8], Data1[9],
                                                Data1[10], Data1[11], Data1[12], Data1[13]);

                Intersection_map.clear();
                ListClusters.clear();
                Percolation_cluster.clear();
                cout << "identifying intersections with truncated fractures" << endl;
                cuDFNsys::IdentifyIntersection<_DataType_> identifyInters2{Frac_verts_host.size(),
                                                                           Frac_verts_device_ptr, true,
                                                                           Intersection_map};
                cout << "identifying cluster with truncated fractures" << endl;
                cuDFNsys::Graph<_DataType_> G2{(size_t)DSIZE, Intersection_map};
                G2.UseDFS(ListClusters);
                cuDFNsys::IdentifyPercolationCluster<_DataType_> IdentiClu2{ListClusters,
                                                                            Frac_verts_host, perco_dir,
                                                                            Percolation_cluster};
                cout << "DFN II finished" << endl;
                cuDFNsys::MatlabPlotDFN<_DataType_> As2{"DFN_II.h5", "DFN_II.m",
                                                        Frac_verts_host, Intersection_map, ListClusters,
                                                        Percolation_cluster, true, true, true, true,
                                                        L, perco_dir, true, "DFN_II", DomainDimensionRatio};
                //return 0;
                Frac_verts_device.clear();
                Frac_verts_device.shrink_to_fit();

                cuDFNsys::OutputObjectData<_DataType_> lk_out;

                if (Percolation_cluster.size() > 0)
                {
                    lk_out.OutputFractures("Fractures.h5", Frac_verts_host, L, DomainDimensionRatio);

                    std::vector<size_t> Fracs_percol;

                    cuDFNsys::GetAllPercolatingFractures GetPer{Percolation_cluster,
                                                                ListClusters,
                                                                Fracs_percol};
                    std::vector<pair<int, int>> IntersectionPair_percol;
                    int NUMprior = Fracs_percol.size();

                    bool ifRemoveDeadEnds = (IfRemoveDeadEnd == 0 ? false : true);
                    cout << "ifRemoveDeadEnds: " << ifRemoveDeadEnds << endl;
                    cuDFNsys::RemoveDeadEndFrac<_DataType_> RDEF{Fracs_percol,
                                                                 IntersectionPair_percol,
                                                                 (size_t)perco_dir,
                                                                 Frac_verts_host,
                                                                 Intersection_map,
                                                                 ifRemoveDeadEnds};

                    cout << "meshing ..." << endl;

                    cuDFNsys::Mesh<_DataType_> mesh;

                    lk_out.OutputFractures("FracturesII.h5", Frac_verts_host, L, DomainDimensionRatio);

                    if (!If_percolative)
                    {
                        cuDFNsys::Mesh<_DataType_> mesh2{Frac_verts_host, IntersectionPair_percol,
                                                         &Fracs_percol, minGridSize, maxGridSize, perco_dir, L,
                                                         DomainDimensionRatio};
                        lk_out.OutputMesh("mesh.h5", mesh2, Fracs_percol);
                        mesh = mesh2;
                        h5g.NewFile("FlowProperties.h5");
                    }
                    else
                        lk.InputMesh("mesh.h5", mesh, &Fracs_percol);

                    double mean_grid_area = mesh.MatlabPlot("DFN_mesh_.h5",
                                                            "DFN_mesh_.m",
                                                            Frac_verts_host,
                                                            L, true, true, true,
                                                            "DFN_mesh_", DomainDimensionRatio);

                    cout << "The mean area of all elements is " << mean_grid_area << endl;

                    cout << "MHFEM ing ..." << endl;

                    cuDFNsys::MHFEM<_DataType_> fem;
                    if (!If_percolative)
                    {
                        cuDFNsys::MHFEM<_DataType_> fem2{mesh, Frac_verts_host, P_in, P_out,
                                                         perco_dir, L, DomainDimensionRatio};
                        lk_out.OutputMHFEM("mhfem.h5", fem2);
                        fem = fem2;
                    }
                    else
                        lk.InputMHFEM("mhfem.h5", fem);

                    cout << "Fluxes: " << fem.QIn << ", ";
                    cout << fem.QOut << ", Permeability: ";
                    cout << fem.Permeability << endl;
                    if (fem.QError > 1 || isnan(fem.Permeability) == 1)
                    {
                        throw cuDFNsys::ExceptionsIgnore("Found large error or isnan, the error: " + std::to_string(fem.QError) + ", the permeability: " + std::to_string(fem.Permeability) + "\n");
                    }
                    double2 TGH = fem.MatlabPlot("MHFEM_.h5",
                                                 "MHFEM_.m",
                                                 Frac_verts_host, mesh, L, true, "MHFEM_", DomainDimensionRatio);

                    double meanV = TGH.x;
                    double maxV = TGH.y;

                    try
                    {
                        h5g.AddDataset<double>("FlowProperties.h5", "N", "MeanV", &meanV, make_uint2(1, 1));
                        h5g.AddDataset<double>("FlowProperties.h5", "N", "maxV", &maxV, make_uint2(1, 1));
                        vector<string> datasetname = {
                            "P33_total",
                            "P33_connected_z",
                            "Ratio_of_P33_z",
                            "P33_largest_cluster",
                            "P32_total",
                            "P32_connected_z",
                            "Ratio_of_P32_z",
                            "P32_largest_cluster",
                            "P30",
                            "P30_connected_z",
                            "Ratio_of_P30_z",
                            "P30_largest_cluster",
                            "Percolation_probability_z",
                            "n_I",
                            "Permeability_z",
                            "Q_error_z"};
                        vector<double *> data_input = {&Data1[0], &Data1[1], &Data1[2], &Data1[3],
                                                       &Data1[4], &Data1[5], &Data1[6], &Data1[7], &Data1[8], &Data1[9],
                                                       &Data1[10], &Data1[11], &Data1[12], &Data1[13], &fem.Permeability, &fem.QError};
                        vector<uint2> dim_ss(data_input.size(), make_uint2(1, 1));
                        h5g.AddDatasetsWithOneGroup<double>("FlowProperties.h5", "ConnectivityPermeability", datasetname, data_input, dim_ss);
                    }
                    catch (H5::Exception &e)
                    {
                        //cout
                        //e.printError();
                        //exit(0);
                        // do nothing
                    }

                    cout << "The maximum velocity of all elements is " << maxV << endl;
                    cout << "The mean velocity of all elements is " << meanV << endl;

                    double meanTime = pow(mean_grid_area, 0.5) / maxV;

                    DeltaT = meanTime / Factor_mean_time_in_grid;
                    cout << "\nThe delta T is set to be " << ("\033[1;33m") << DeltaT << ("\033[0m") << "\n\n";

                    DiffusionLocal = LengthScale_Over_Pe * meanV;
                    cout << "\nThe DiffusionLocal is set to be " << ("\033[1;33m") << DiffusionLocal << ("\033[0m") << "\n\n";
                    //-----------------

                    string FractureFileName_r = "ParticlePositionResult/DispersionInfo.h5";

                    std::ifstream fileeqr(FractureFileName_r);
                    bool psd = fileeqr.good();

                    if (psd)
                    {
                        cuDFNsys::HDF5API hg6;
                        vector<double> Aqs = hg6.ReadDataset<double>(FractureFileName_r,
                                                                     "N", "Delta_T");
                        cout << "\nThe delta T is set to be " << ("\033[1;33m") << Aqs[0] << ("\033[0m") << "\n\n";

                        Aqs = hg6.ReadDataset<double>(FractureFileName_r,
                                                      "N", "Dispersion_local");
                        cout << "\nThe DiffusionLocal is set to be " << ("\033[1;33m") << Aqs[0] << ("\033[0m") << "\n\n";
                    }

                    string Filename_FracturesForParticle = "FracturesForParticle.h5";

                    std::ifstream file(Filename_FracturesForParticle);
                    bool pwqs = file.good();

                    if (!pwqs)
                    {
                        cout << "Writting " << Filename_FracturesForParticle << endl;
                        lk_out.OutputFractures(Filename_FracturesForParticle, Frac_verts_host, L);
                    }
                    cout << "Particle transport ing ......\n";

                    cuDFNsys::ParticleTransport<_DataType_> p{NumTimeSteps_Dispersion, // number of time step
                                                              Frac_verts_host, mesh, fem, (uint)perco_dir,
                                                              -0.5f * L * (&DomainDimensionRatio.x)[perco_dir],
                                                              NumParticlesRandomWalk, // num of particle
                                                              DeltaT,                 // delta_T_ii
                                                              DiffusionLocal,
                                                              "Particle_tracking",
                                                              "Flux-weighted",
                                                              recordMode,
                                                              false, 1, false, ControlPlaneSpacing, IfoutputMsd};
                    p.MatlabPlot("MHFEM_.h5", "ParticlesDFNMatlab.m", mesh, fem, L, DomainDimensionRatio, true, "ParticlesDFN");
                }
                else
                {
                    system("rm -rf *.h5 ParticlePositionResult");
                    j = 0;
                }
            }
            catch (cuDFNsys::ExceptionsIgnore &e)
            {
                cout << "cuDFNsys::ExceptionsIgnore\n";
                cout << e.what() << endl;
                cout << path2 << endl;
                system("rm -rf *.h5 ParticlePositionResult");
                j = 0;
            }
            catch (cuDFNsys::ExceptionsPause &e)
            {
                cout << "cuDFNsys::ExceptionsPause\n";
                cout << e.what() << endl;
                cout << path2 << endl;
                system("rm -rf *.h5 ParticlePositionResult");
                j = 0;
            }
            catch (H5::Exception &e)
            {
                cout << "H5::Exception\n";
                //e.printError();
                cout << path2 << endl;
                system("rm -rf *.h5 ParticlePositionResult");
                j = 0;
            }
            catch (...)
            {
                cout << "Unknown exceptions!\n";
                cout << path2 << endl;
                system("rm -rf *.h5 ParticlePositionResult");
                j = 0;
            }
        }

        chdir(curPath.c_str());
    }

    //-------------------------------------------

    return 0;
};

int GetRowsOfDataset(string Filename, string Datasetname)
{
    H5File file(Filename, H5F_ACC_RDONLY);
    DataSet dataset;

    dataset = file.openDataSet(Datasetname);

    DataSpace filespace = dataset.getSpace();
    int rank = filespace.getSimpleExtentNdims();

    hsize_t dims[rank];
    rank = filespace.getSimpleExtentDims(dims);

    file.close();

    return int(dims[1]);
};