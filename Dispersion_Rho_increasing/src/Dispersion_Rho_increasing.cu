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
// DESCRIPTION: Dispersion in a DFN with a specific percolation parameter value
// AUTHOR:      Tingchang YIN
// DATE:        24/03/2023
// ====================================================

#include "cuDFNsys.cuh"
#include <fstream>
#include <iostream>
#include <limits.h>
#include <unistd.h>
#ifdef USE_DOUBLES
typedef double _DataType_;
#else
typedef float _DataType_;
#endif

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
    uint LoopTimes = atoi(argv[1]);
    uint InitDensity = atoi(argv[2]);
    uint DensityIncreament = atoi(argv[3]);

    //------------ other inputs
    _DataType_ L = atof(argv[4]);
    _DataType_ kappa_ = atof(argv[5]),
               beta_ = atof(argv[6]),
               gamma_ = atof(argv[7]);
    int size_frac_mode = atoi(argv[8]); // mode of fracture size distributions
    cuDFNsys::Vector4<_DataType_> ParaSizeDistri =
        cuDFNsys::MakeVector4((_DataType_)atof(argv[9]),
                              (_DataType_)atof(argv[10]),
                              (_DataType_)atof(argv[11]),
                              (_DataType_)atof(argv[12]));
    double3 DomainDimensionRatio = make_double3(1, 1, atof(argv[13]));
    int IfRemoveDeadEnd = atoi(argv[14]);
    _DataType_ minGridSize = atof(argv[15]);
    _DataType_ maxGridSize = atof(argv[16]);

    int NumTimeSteps_Dispersion = atoi(argv[17]);
    int NumParticlesRandomWalk = atoi(argv[18]);
    _DataType_ DeltaT = 0;
    _DataType_ Factor_mean_time_in_grid = atof(argv[19]);
    // the mean time (a characteristic grid length over the mean velocity (m/s)) for a random walk to cross a characteristic grid length
    // but this mean time was reduced, i.e., dividing by a factor (> 1)
    // then the mean time is DeltaT
    _DataType_ DiffusionLocal = 0;
    _DataType_ LengthScale_Over_Pe = 0;
    _DataType_ LengthScale = atof(argv[20]);
    _DataType_ Pe = atof(argv[21]);
    _DataType_ ControlPlaneSpacing = atof(argv[22]);
    bool IfoutputMsd = atoi(argv[23]) == 0 ? false : true;
    bool IfoutputParticleInfoAllsteps = atoi(argv[24]) == 0 ? false : true;
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

    for (uint i = 1; i <= LoopTimes; ++i)
    {
        string path2 = "DFN_" + cuDFNsys::ToStringWithWidth(i, 3);
        string command1 = "mkdir -p " + path2;
        system(command1.c_str());

        string command2 = curPath + "/" + path2;
        chdir(command2.c_str());

        uint DSIZE = InitDensity + (i - 1) * DensityIncreament;

        for (uint j = 0;; j++)
        {
            //-----------if a mesh h5 exists ----
            // then the DFN is percolative
            bool If_percolative = false;
            std::ifstream fileer("mesh.h5");
            If_percolative = fileer.good();

            //----------if a mesh exists
            //---------- the dispersion simultion may have been finished
            //------ let's check it
            if(If_percolative)
            {
                std::ifstream FileYe("ParticlePositionResult/ParticlePositionLastStep.h5"); 
                bool DFSW = FileYe.good(); 

                if (DFSW) // that means the simulation has been conducted
                {

                }
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
                cuDFNsys::Fractures<_DataType_><<<DSIZE / 256 + 1, 256>>>(Frac_verts_device_ptr,
                                                                          (unsigned long)t + (unsigned long)ceil(abs(Ter(0, 0)) * ((unsigned long)t * 1.0)),
                                                                          DSIZE, L,
                                                                          0,
                                                                          ParaSizeDistri,
                                                                          kappa_, // kappa
                                                                          beta_,  // beta
                                                                          gamma_, // gamma
                                                                          DomainDimensionRatio);
            else
            {
                Frac_verts_host.resize(0);
                lk.InputFractures("Fractures.h5", Frac_verts_host, L, DomainDimensionRatio);
                Frac_verts_host.shrink_to_fit();
            }
            cudaDeviceSynchronize();
        }

        chdir(curPath.c_str());
    }

    //-------------------------------------------

    return 0;
};
