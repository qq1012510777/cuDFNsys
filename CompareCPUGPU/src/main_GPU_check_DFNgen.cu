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
// NAME:        CompareCPUGPU
// DESCRIPTION: CompareCPUGPU
// AUTHOR:      Tingchang YIN
// DATE:        40/11/2022
// ====================================================

#include "cuDFNsys.cuh"
#include <bits/stdc++.h>
#include <iostream>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#ifdef USE_DOUBLES
typedef double _DataType_;
#else
typedef float _DataType_;
#endif

int main(int argc, char *argv[])
{

    try
    {
        int dev = 0;
        GPUErrCheck(cudaSetDevice(dev));

        cuDFNsys::Warmup<<<1024 / 256 + 1, 256 /*  1, 2*/>>>();
        cudaDeviceSynchronize();

        _DataType_ L = 30;
        //int perco_dir = 2;
        cuDFNsys::Vector4<_DataType_> sizedistribution_para = cuDFNsys::MakeVector4(1.5, 1., 15., 0.);

        uint StartLoopStep = atoi(argv[1]);
        uint EndLoopStep = atoi(argv[2]);
        uint MonteCarloTimes = atoi(argv[3]);
        int DensityIncrement = atoi(argv[4]);

        string filename = "GPU_check_DFNgen";
        cuDFNsys::HDF5API hd5;

        string hdfilename = filename + "/" + "DFN_gen_countTime.h5";
        uint2 Dims = make_uint2(1, MonteCarloTimes);

        for (uint i = StartLoopStep; i <= EndLoopStep; ++i)
        {
            cout << "step " << i << endl;
            if (i == 1)
            {
                int us = remove(filename.c_str());

                us = mkdir(filename.c_str(), 0777);

                if (us == -1)
                    throw cuDFNsys::ExceptionsPause("Cannot create file " + filename);

                hd5.NewFile(hdfilename);

                hd5.AddDataset(hdfilename, "N",
                               "DensityIncrement", &DensityIncrement, make_uint2(1, 1));
                hd5.AddDataset(hdfilename, "N",
                               "EndLoopStep", &EndLoopStep, make_uint2(1, 1));
                hd5.AddDataset(hdfilename, "N",
                               "MonteCarloTimes", &MonteCarloTimes, make_uint2(1, 1));
            }

            vector<double> DFNgenTime(MonteCarloTimes);
            vector<double> n_I(MonteCarloTimes);

            double istart_all = cuDFNsys::CPUSecond();
            for (uint j = 0; j < MonteCarloTimes; ++j)
            {

                time_t t;
                time(&t);
                t = t * (j + 1);

                //cout << "\t\t" << (unsigned long)t << endl; 

                int DSIZE = i * DensityIncrement;

                thrust::host_vector<cuDFNsys::Fracture<_DataType_>> Frac_verts_host(DSIZE);
                thrust::device_vector<cuDFNsys::Fracture<_DataType_>> Frac_verts_device(DSIZE);
                cuDFNsys::Fracture<_DataType_> *Frac_verts_device_ptr;
                Frac_verts_device_ptr = thrust::raw_pointer_cast(Frac_verts_device.data());

                double istart = cuDFNsys::CPUSecond();
                cuDFNsys::Fractures<_DataType_><<<DSIZE / 256 + 1, 256>>>(Frac_verts_device_ptr,
                                                                          (unsigned long)t,
                                                                          DSIZE,
                                                                          L,
                                                                          (uint)0,               // 0 = power law; 1 = lognormal; 2 = uniform; 3 = monosize
                                                                          sizedistribution_para, // when mode = 1, ;
                                                                          (_DataType_)0.0,
                                                                          (_DataType_)0.25);

                std::map<pair<size_t, size_t>, pair<cuDFNsys::Vector3<_DataType_>, cuDFNsys::Vector3<_DataType_>>> Intersection_map;
                cudaDeviceSynchronize();

                Frac_verts_host = Frac_verts_device;

                cuDFNsys::IdentifyIntersection<_DataType_> identifyInters{Frac_verts_host.size(),
                                                                          Frac_verts_device_ptr,
                                                                          false,
                                                                          Intersection_map};
                double elapse_ = cuDFNsys::CPUSecond() - istart;

                DFNgenTime[j] = elapse_;
                n_I[j] = 1.0f * Intersection_map.size() / (Frac_verts_host.size() * 1.0f);
                // cout << "\t\t" << n_I[j] << endl;
            }

            double ielapse_all = cuDFNsys::CPUSecond() - istart_all;

            hd5.AddDataset(hdfilename, "Step_" + cuDFNsys::ToStringWithWidth(i, 5),
                           "DFNgenTime", DFNgenTime.data(), Dims);
            hd5.AddDataset(hdfilename, "Step_" + cuDFNsys::ToStringWithWidth(i, 5),
                           "n_I", n_I.data(), Dims);

            hd5.AddDataset(hdfilename, "Step_" + cuDFNsys::ToStringWithWidth(i, 5),
                           "SumOfTime", &ielapse_all, make_uint2(1, 1));
        }
    }
    catch (cuDFNsys::ExceptionsIgnore &e)
    {
        cout << e.what() << endl;
    }
    catch (cuDFNsys::ExceptionsPause &e)
    {
        cout << e.what() << endl;
        return 0;
    }
    catch (...)
    {
        throw;
    };
    return 0;
};