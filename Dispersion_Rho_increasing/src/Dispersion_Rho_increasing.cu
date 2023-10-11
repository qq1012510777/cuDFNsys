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
#include <numeric>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>
#ifdef USE_DOUBLES
typedef double _DataType_;
#else
typedef float _DataType_;
#endif

int GetRowsOfDataset(string Filename, string Datasetname);
string GetArgvWithoutComment(const string &A, const char &B);

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
    char split_char = '#';

    uint InitLoop = atoi(GetArgvWithoutComment(argv[1], split_char).c_str());
    uint FinalLoop = atoi(GetArgvWithoutComment(argv[2], split_char).c_str());
    uint InitDensity = atoi(GetArgvWithoutComment(argv[3], split_char).c_str());
    uint DensityIncreament = atoi(GetArgvWithoutComment(argv[4], split_char).c_str());
    uint MaxTranLoopTimes = atoi(GetArgvWithoutComment(argv[5], split_char).c_str());

    uint LoopTimes = FinalLoop - InitLoop + 1;

    //------------ other inputs
    _DataType_ L = atof(GetArgvWithoutComment(argv[6], split_char).c_str());
    _DataType_ kappa_ = atof(GetArgvWithoutComment(argv[7], split_char).c_str()),
               beta_ = atof(GetArgvWithoutComment(argv[8], split_char).c_str()),
               gamma_ = atof(GetArgvWithoutComment(argv[9], split_char).c_str());
    int size_frac_mode = atoi(GetArgvWithoutComment(argv[10], split_char).c_str()); // mode of fracture size distributions
    cuDFNsys::Vector4<_DataType_> ParaSizeDistri =
        cuDFNsys::MakeVector4((_DataType_)atof(GetArgvWithoutComment(argv[11], split_char).c_str()),
                              (_DataType_)atof(GetArgvWithoutComment(argv[12], split_char).c_str()),
                              (_DataType_)atof(GetArgvWithoutComment(argv[13], split_char).c_str()),
                              (_DataType_)atof(GetArgvWithoutComment(argv[14], split_char).c_str()));
    double3 DomainDimensionRatio = make_double3(1, 1, atof(GetArgvWithoutComment(argv[15], split_char).c_str()));
    int IfRemoveDeadEnd = atoi(GetArgvWithoutComment(argv[16], split_char).c_str());
    _DataType_ minGridSize = atof(GetArgvWithoutComment(argv[17], split_char).c_str());
    _DataType_ maxGridSize = atof(GetArgvWithoutComment(argv[18], split_char).c_str());

    int NumTimeSteps_Dispersion = atoi(GetArgvWithoutComment(argv[19], split_char).c_str());
    int NumParticlesRandomWalk = atoi(GetArgvWithoutComment(argv[20], split_char).c_str());
    _DataType_ DeltaT = 0;
    _DataType_ Factor_mean_time_in_grid = atof(GetArgvWithoutComment(argv[21], split_char).c_str());
    // the mean time (a characteristic grid length over the mean velocity (m/s)) for a random walk to cross a characteristic grid length
    // but this mean time was reduced, i.e., dividing by a factor (> 1)
    // then the mean time is DeltaT
    _DataType_ DiffusionLocal = 0;
    _DataType_ LengthScale_Over_Pe = 0;
    _DataType_ LengthScale = atof(GetArgvWithoutComment(argv[22], split_char).c_str());
    _DataType_ Pe = atof(GetArgvWithoutComment(argv[23], split_char).c_str());
    _DataType_ ControlPlaneSpacing = atof(GetArgvWithoutComment(argv[24], split_char).c_str());
    bool IfoutputMsd = atoi(GetArgvWithoutComment(argv[25], split_char).c_str()) == 0 ? false : true;
    bool IfoutputParticleInfoAllsteps = atoi(GetArgvWithoutComment(argv[26], split_char).c_str()) == 0 ? false : true;
    int ThresholdToStop = atoi(GetArgvWithoutComment(argv[27], split_char).c_str());
    int ThresholdForMaximumLeftParticles = atoi(GetArgvWithoutComment(argv[28], split_char).c_str());

    string injectionMode = GetArgvWithoutComment(argv[29], split_char).c_str(); // Flux-weighted or Resident
    bool IfInjectAt_Center = (atoi(GetArgvWithoutComment(argv[30], split_char).c_str()) == 0 ? false : true);
    _DataType_ InjectionPlane = atof(GetArgvWithoutComment(argv[31], split_char).c_str());
    bool If_ReRun_ = atoi(GetArgvWithoutComment(argv[32], split_char).c_str()) == 0 ? false : true;
    string LogFile = GetArgvWithoutComment(argv[33], split_char).c_str();
    bool IfComplexeMixing = atoi(GetArgvWithoutComment(argv[34], split_char).c_str()) == 0 ? false : true;

    bool IfStopAtFirstArrival = false;
    if (GetArgvWithoutComment(argv[35], split_char).c_str() != NULL)
        IfStopAtFirstArrival = (atoi(GetArgvWithoutComment(argv[35], split_char).c_str()) == 0 ? false : true);
    else
        IfStopAtFirstArrival = false;

    bool IfUseMeanV_or_DarcianQ_toCalculatedPe = true;
    if (GetArgvWithoutComment(argv[36], split_char).c_str() != NULL)
        IfUseMeanV_or_DarcianQ_toCalculatedPe = (atoi(GetArgvWithoutComment(argv[36], split_char).c_str()) == 0 ? false : true);
    else
        IfUseMeanV_or_DarcianQ_toCalculatedPe = true;

    uint ReRunMaxTime = 0;
    if (GetArgvWithoutComment(argv[37], split_char).c_str() != NULL)
        ReRunMaxTime = atoi(GetArgvWithoutComment(argv[37], split_char).c_str());
    else
        ReRunMaxTime = 5;

    string recordMode = IfoutputParticleInfoAllsteps == false ? "FPTCurve" : "OutputAll";
    _DataType_ P_in = L, P_out = 0;

    std::ofstream oss("./InputInfo.info", ios::out);

    oss << "L: " << L << endl;
    oss << "Kappa: " << kappa_ << endl;
    oss << "Beta: " << beta_ << endl;
    oss << "Gamma: " << gamma_ << endl;
    oss << "Mode of fracture size distributions: " << size_frac_mode << endl;
    oss << "Parameters of the size distribution: " << ParaSizeDistri.x << ", " << ParaSizeDistri.y << ", " << ParaSizeDistri.z << ", " << ParaSizeDistri.w << endl;
    oss << "Domain's dimension ratio: " << DomainDimensionRatio.x << ", " << DomainDimensionRatio.y << ", " << DomainDimensionRatio.z << endl;
    oss << "If remove the dead ends: " << (IfRemoveDeadEnd == 0 ? "false" : "true") << endl;
    oss << "Min grid size: " << minGridSize << endl;
    oss << "Max grid size: " << maxGridSize << endl;
    oss << "Hydraulic head at the inlet and outlet: " << P_in << ", " << P_out << endl;
    oss << "Number of time steps for random walks: " << NumTimeSteps_Dispersion << endl;
    oss << "Number of particles: " << NumParticlesRandomWalk << endl;
    oss << "Factor_mean_time_in_grid: " << Factor_mean_time_in_grid << endl;
    oss << "LengthScale: " << LengthScale << endl;
    oss << "Pe: " << Pe << endl;
    oss << "The spacing of control planes: " << ControlPlaneSpacing << endl;
    oss << "IfoutputMsd: " << (IfoutputMsd == true ? "true" : "false") << endl;
    oss << "IfoutputParticleInfoAllsteps: " << (IfoutputParticleInfoAllsteps == false ? "FPTCurve" : "OutputAll") << endl;
    oss << "injectionMode: " << injectionMode << endl;
    oss << "InjectionPlane: " << (IfInjectAt_Center ? InjectionPlane : 0.5 * L) << endl;
    oss << "If_ReRun_: " << (If_ReRun_ ? "true" : "false") << endl;
    oss << "IfComplexeMixing: " << (IfComplexeMixing ? "Outlet-flux-weighted" : "Equal-probability") << endl;
    oss << "IfUseMeanV_or_DarcianQ_toCalculatedPe: " << (IfUseMeanV_or_DarcianQ_toCalculatedPe ? "MeanV" : "DarcianQ") << endl;
    oss << "ReRunMaxTime: " << ReRunMaxTime << endl;
    oss << "IfStopAtFirstArrival: " << (IfStopAtFirstArrival ? "true" : "false") << endl;
    oss.close();

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

        if (i == InitLoop)
            system("echo \" \" > ../recordTime_and_Error.log");

        string Start_i = "loop " + std::to_string(i) + ", NF = " + std::to_string(DSIZE) + " started\n";
        system("echo $(date +%d-%m-%y---%T) >> ../recordTime_and_Error.log");
        Start_i = "echo \"" + Start_i + "\" >> ../recordTime_and_Error.log";
        system(Start_i.c_str());

        uint countReRunTime = 0; // ReRunMaxTime
        for (uint j = 0; j < MaxTranLoopTimes; j++)
        {

            try
            {
                string GGF = "echo \" \" > ../" + LogFile;

                system(GGF.c_str());

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
                    std::ifstream FileYe("ParticlePositionResult/ParticlePosition_WhichStepDoesTheParticleReached.h5");
                    bool DFSW = FileYe.good();

                    if (DFSW) // that means the simulation has been conducted
                    {
                        std::vector<double> DF = h5g.ReadDataset<double>("ParticlePositionResult/ParticlePosition_WhichStepDoesTheParticleReached.h5", "N",
                                                                         "WhichStepDoesTheParticleReached");
                        //----------------------
                        if (IfStopAtFirstArrival)
                        {
                            double resusdslt = accumulate(DF.begin(), DF.end(), 0);
                            if (resusdslt != -1.0 * DF.size())
                            {
                                system("echo $(date +%d-%m-%y---%T) >> ../recordTime_and_Error.log");
                                system("echo \"Finished\n\n\" >> ../recordTime_and_Error.log");
                                break; // get enough particles arrived  }
                            }
                        }

                        //-----------------------------------
                        std::vector<uint> NumPaLeft = h5g.ReadDataset<uint>("ParticlePositionResult/DispersionInfo.h5",
                                                                            "N", "NumParticlesLeftFromInlet");
                        if (NumPaLeft[0] > ThresholdForMaximumLeftParticles)
                        {
                            system("echo $(date +%d-%m-%y---%T) >> ../recordTime_and_Error.log");
                            system("echo \"Too many particles left from the inlet\n\n\" >> ../recordTime_and_Error.log");
                            throw cuDFNsys::ExceptionsIgnore("Too many particles left from the inlet!\n");
                        }

                        //-----------------------------------
                        if (!IfoutputParticleInfoAllsteps)
                        {
                            std::vector<uint> NUMsTEPS = h5g.ReadDataset<uint>("ParticlePositionResult/DispersionInfo.h5",
                                                                               "N", "NumOfSteps");

                            int rows_ = GetRowsOfDataset("ParticlePositionResult/ParticlePositionLastStep.h5",
                                                         "Step_" + cuDFNsys::ToStringWithWidth<int>(NUMsTEPS[0], 10));
                            if (rows_ <= ThresholdToStop)
                            {
                                system("echo $(date +%d-%m-%y---%T) >> ../recordTime_and_Error.log");
                                system("echo \"Finished\n\n\" >> ../recordTime_and_Error.log");
                                break; // get enough particles arrived
                            }
                        }
                        else
                        {
                            int NumParticles = DF.size() / 2;
                            // NumPaLeft

                            int NumParArrival = NumParticles - count(DF.begin(), DF.begin() + NumParticles, -1);

                            if ((NumParticles - NumParArrival - NumPaLeft[0]) <= ThresholdToStop)
                            {
                                system("echo $(date +%d-%m-%y---%T) >> ../recordTime_and_Error.log");
                                system("echo \"Finished\n\n\" >> ../recordTime_and_Error.log");
                                break; // get enough particles arrived
                            }
                        }
                    }

                    //----------amend MSD data
                    system("python3 ../scripts/DelDataSet.py");
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

                    if (!IfUseMeanV_or_DarcianQ_toCalculatedPe)
                        DiffusionLocal = LengthScale_Over_Pe * fem.Permeability; //meanV; // change to Darcian q
                    else
                        DiffusionLocal = LengthScale_Over_Pe * meanV;

                    //-----------------------------------advction time to cross a grid
                    _DataType_ DeltaT_Advection = pow(mean_grid_area, 0.5) / maxV;
                    DeltaT_Advection = DeltaT_Advection / Factor_mean_time_in_grid;
                    //-----------------------------------diffusion time to cross a grid
                    _DataType_ DeltaT_Diffusion = pow(mean_grid_area, 0.5) / (pow(2 * DiffusionLocal, 0.5) * 4);
                    DeltaT_Diffusion = DeltaT_Diffusion * DeltaT_Diffusion / Factor_mean_time_in_grid;
                    cout << "DeltaT_Advection: " << DeltaT_Advection << ", DeltaT_Diffusion: " << DeltaT_Diffusion << endl;
                    if (DeltaT_Advection < DeltaT_Diffusion)
                    {
                        cout << "\n** ADVEction Time is smaller **\n\n";
                        DeltaT = DeltaT_Advection;
                    }
                    else
                    {
                        cout << "\n** DIFFusion Time is smaller **\n\n";
                        DeltaT = DeltaT_Diffusion;
                    }

                    // if (IfComplexeMixing) // advection dominated
                    // {
                    //     cout << "** advection dominated **\n";
                    //     double meanTime = pow(mean_grid_area, 0.5) / maxV;
                    //     DeltaT = meanTime / Factor_mean_time_in_grid;
                    // }
                    // else
                    // {
                    //     cout << "** diffusion dominated **\n";
                    //     DeltaT = pow(mean_grid_area, 0.5) / (pow(2 * DiffusionLocal, 0.5) * Factor_mean_time_in_grid);
                    //     DeltaT = DeltaT * DeltaT;
                    // }

                    cout << "The maximum velocity of all elements is " << maxV << endl;
                    cout << "The mean velocity of all elements is " << meanV << endl;

                    cout << "\nThe delta T is set to be " << ("\033[1;33m") << DeltaT << ("\033[0m") << "\n\n";

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
                                                              injectionMode,
                                                              recordMode,
                                                              false, 1, false, ControlPlaneSpacing, IfoutputMsd,
                                                              IfInjectAt_Center,
                                                              InjectionPlane, IfComplexeMixing};
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

                system("echo $(date +%d-%m-%y---%T) >> ../recordTime_and_Error.log");
                string ERDS = "echo \"\t" + e.what() + "\n\n\" >> ../recordTime_and_Error.log";
                system(ERDS.c_str());
            }
            catch (cuDFNsys::ExceptionsPause &e)
            {
                cout << "cuDFNsys::ExceptionsPause\n";
                cout << e.what() << endl;
                cout << path2 << endl;
                if (!If_ReRun_)
                {
                    system("rm -rf *.h5 ParticlePositionResult");
                    j = 0;
                }
                else
                { // re-run
                    string vccy = "result=$(grep \"Loop times is too large!\" ../" + LogFile + ")";
                    int status = system(vccy.c_str());

                    if (status == 0)
                    {
                        countReRunTime++;

                        if (countReRunTime <= ReRunMaxTime)
                        {
                            system("echo $(date +%d-%m-%y---%T) >> ../recordTime_and_Error.log");
                            string ERDS = "echo \"ReRunTheSameDFN_Transport\" >> ../recordTime_and_Error.log";
                            system(ERDS.c_str());
                        }
                        else
                        {
                            string ERDS = "echo \"ReRunTimes reached the maximum\" >> ../recordTime_and_Error.log";
                            system(ERDS.c_str());

                            system("rm -rf *.h5 ParticlePositionResult");
                            j = 0;
                        }
                        //exit(0);
                    }
                    else
                    {
                        system("rm -rf *.h5 ParticlePositionResult");
                        j = 0;
                        //exit(0);
                    }
                }

                system("echo $(date +%d-%m-%y---%T) >> ../recordTime_and_Error.log");
                string ERDS = "echo \"\t" + e.what() + "\n\n\" >> ../recordTime_and_Error.log";
                system(ERDS.c_str());
            }
            catch (H5::Exception &e)
            {
                cout << "H5::Exception\n";
                //e.printError();
                cout << path2 << endl;
                system("rm -rf *.h5 ParticlePositionResult");
                j = 0;

                system("echo $(date +%d-%m-%y---%T) >> ../recordTime_and_Error.log");
                system("echo \"H5::Exception\n\n\" >> ../recordTime_and_Error.log");
            }
            catch (...)
            {
                cout << "Unknown exceptions!\n";
                cout << path2 << endl;
                system("rm -rf *.h5 ParticlePositionResult");
                j = 0;

                system("echo $(date +%d-%m-%y---%T) >> ../recordTime_and_Error.log");
                system("echo \"Unknown exceptions\n\n\" >> ../recordTime_and_Error.log");
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

string GetArgvWithoutComment(const string &A, const char &B)
{
    std::stringstream test(A);
    std::string segment;
    std::vector<std::string> seglist;

    while (std::getline(test, segment, B))
    {
        seglist.push_back(segment);
    }

    return seglist[0];
};