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
// NAME:        A test case
// DESCRIPTION: Call cuDFNsys functions to do simulation and test.
// AUTHOR:      Tingchang YIN
// DATE:        30/06/2022
// ====================================================

#include "cuDFNsys.cuh"
#include <omp.h>
#include <unistd.h>

#ifdef USE_DOUBLES
typedef double _DataType_;
#else
typedef float _DataType_;
#endif

__global__ void Transform2DTO3DKernel(cuDFNsys::Fracture<_DataType_> *Frac_verts_device_ptr,
                                      _DataType_ *Position3D_dev_ptr,
                                      _DataType_ *temp2Dpos_dev_ptr,
                                      uint *ElementFracTag_cuda_devptr,
                                      uint *EleTag_device_ptr,
                                      uint count)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;

    if (idx > count - 1)
        return;

    uint numParticles = count;
    cuDFNsys::Vector3<_DataType_> Pos = cuDFNsys::MakeVector3(temp2Dpos_dev_ptr[idx],
                                                              temp2Dpos_dev_ptr[idx + numParticles],
                                                              (_DataType_)0.0);
    //printf("temp2Dpos_dev_ptr: %f, %f\n", Pos.x, Pos.y);
    uint EleTag_j = EleTag_device_ptr[idx];

    uint FracTag_j = ElementFracTag_cuda_devptr[EleTag_j - 1] - 1;

    //printf("%d, FracTag_j: %d\n", idx, FracTag_j);

    _DataType_ Rotate2DTo3D[3][3];
    Frac_verts_device_ptr[FracTag_j].RoationMatrix(Rotate2DTo3D, 23);

    Pos = cuDFNsys::ProductSquare3Vector3<_DataType_>(Rotate2DTo3D, Pos);
    Pos = cuDFNsys::MakeVector3(Pos.x + Frac_verts_device_ptr[FracTag_j].Center.x,
                                Pos.y + Frac_verts_device_ptr[FracTag_j].Center.y,
                                Pos.z + Frac_verts_device_ptr[FracTag_j].Center.z);

    Position3D_dev_ptr[idx] = Pos.x;
    Position3D_dev_ptr[idx + numParticles] = Pos.y;
    Position3D_dev_ptr[idx + 2 * numParticles] = Pos.z;
};

int main(int argc, char *argv[])
{

    try
    {
        _DataType_ L;
        //uint DSIZE;
        //uint Nproc = 10;
        //if (argv[1] != NULL)
        //    Nproc = atoi(argv[1]);

        string FracH5 = "FracturesForParticle.h5";

        if (argv[1] == NULL)
            throw cuDFNsys::ExceptionsPause("Please give the name of mesh file (.h5)!");

        string mshfile = argv[1];

        thrust::host_vector<cuDFNsys::Fracture<_DataType_>> Frac_verts_host;

        cuDFNsys::InputObjectData<_DataType_> lk;
        lk.InputFractures(FracH5, Frac_verts_host, L);
        //DSIZE = Frac_verts_host.size();
        thrust::device_vector<cuDFNsys::Fracture<_DataType_>> Frac_verts_device;
        Frac_verts_device = Frac_verts_host;
        cuDFNsys::Fracture<_DataType_> *Frac_verts_device_ptr;
        Frac_verts_device_ptr = thrust::raw_pointer_cast(Frac_verts_device.data());

        cuDFNsys::Warmup<<<256 / 256 + 1, 256 /*  1, 2*/>>>();
        cudaDeviceSynchronize();

        string ParticlePosition = "ParticlePositionResult/ParticlePosition";
        string DispersionInfo = "ParticlePositionResult/DispersionInfo";

        string DispersionInfo_h5 = DispersionInfo + ".h5";

        cuDFNsys::HDF5API h5g;
        vector<uint> Tem_p = h5g.ReadDataset<uint>(DispersionInfo_h5, "N", "NumOfSteps");
        uint ExistingNumsteps = Tem_p[0];

        Tem_p = h5g.ReadDataset<uint>(DispersionInfo_h5, "N", "BlockNOPresent");
        uint BlockNOPresent = Tem_p[0];

        Tem_p = h5g.ReadDataset<uint>(DispersionInfo_h5, "N", "SizeOfDataBlock");
        uint SizeOfDataBlock = Tem_p[0];

        vector<uint> ElementFracTag = h5g.ReadDataset<uint>(mshfile, "N", "element_Frac_Tag");

        thrust::host_vector<uint> ElementFracTag_cuda(ElementFracTag.size());
        for (uint i = 0; i < ElementFracTag.size(); ++i)
            ElementFracTag_cuda[i] = ElementFracTag[i];
        thrust::device_vector<uint> ElementFracTag_cuda_dev = ElementFracTag_cuda;
        uint *ElementFracTag_cuda_devptr = thrust::raw_pointer_cast(ElementFracTag_cuda_dev.data());

        //cout << ElementFracTag << endl;
        for (uint i = 0; i <= ExistingNumsteps; ++i)
        {
            cout << "step " << i << " is outputing ...\n";
            string filename = " ", outfilename = "";
            if (i == 0)
            {

                filename = ParticlePosition + "Init.h5";
                outfilename = ParticlePosition + "Init_3D.h5";
                h5g.NewFile(outfilename);
            }
            else
            {
                filename = ParticlePosition + "Block" +
                           cuDFNsys::ToStringWithWidth((uint)ceil((double)i / ((double)SizeOfDataBlock)), 10) + ".h5";

                outfilename = ParticlePosition + "Block" +
                              cuDFNsys::ToStringWithWidth((uint)ceil((double)i / ((double)SizeOfDataBlock)), 10) + "_3D.h5";
                if (i % SizeOfDataBlock == 1)
                    h5g.NewFile(outfilename);
            }

            vector<_DataType_> temp2Dpos = h5g.ReadDataset<_DataType_>(filename, "N", "Step_" + cuDFNsys::ToStringWithWidth(i, 10));
            _DataType_ *data_s = temp2Dpos.data();
            thrust::host_vector<_DataType_> temp2DposCUDA(data_s, data_s + temp2Dpos.size());
            thrust::device_vector<_DataType_> temp2DposCUDA_dev = temp2DposCUDA;
            _DataType_ *temp2DposCUDA_dev_ptr;
            temp2DposCUDA_dev_ptr = thrust::raw_pointer_cast(temp2DposCUDA_dev.data());
            uint numParticles = temp2Dpos.size() / 2;
            //cout << "numParticles: " << numParticles << endl;

            vector<uint> EleTag = h5g.ReadDataset<uint>(filename, "N", "ParticleIDAndElementTag_" + cuDFNsys::ToStringWithWidth(i, 10));
            uint *data_y = EleTag.data();
            thrust::host_vector<uint> EleTag_host(data_y + EleTag.size() / 2, data_y + EleTag.size());
            /// for (uint k = 0; k < EleTag_host.size(); ++k)
            /// cout << EleTag_host[k] << (k < EleTag_host.size() - 1 ? ", " : "\n");
            thrust::device_vector<uint> EleTag_device = EleTag_host;
            uint *EleTag_device_ptr = thrust::raw_pointer_cast(EleTag_device.data());

            //_DataType_ *Position3D = new _DataType_[numParticles * 3];

            thrust::host_vector<_DataType_> Position3D(numParticles * 3);
            thrust::device_vector<_DataType_> Position3D_dev = Position3D;
            _DataType_ *Position3D_dev_ptr = thrust::raw_pointer_cast(Position3D_dev.data());

            Transform2DTO3DKernel<<<numParticles / 256 + 1, 256 /*  1, 2*/>>>(Frac_verts_device_ptr,
                                                                              Position3D_dev_ptr,
                                                                              temp2DposCUDA_dev_ptr,
                                                                              ElementFracTag_cuda_devptr,
                                                                              EleTag_device_ptr,
                                                                              numParticles);
            cudaDeviceSynchronize();
            Position3D = Position3D_dev;
            //#pragma omp parallel for schedule(static) num_threads(Nproc)
            //for (uint j = 0; j < numParticles; ++j)
            //{
            //                cuDFNsys::Vector3<_DataType_> Pos = cuDFNsys::MakeVector3(temp2Dpos[j], temp2Dpos[j + numParticles], (_DataType_)0.0);
            //                //cout << "2D:  " << temp2Dpos[j] << ", " << temp2Dpos[j + numParticles] << endl;
            //
            //                uint EleTag_j = EleTag[numParticles + j];
            //
            //                uint FracTag_j = ElementFracTag(EleTag_j - 1, 0) - 1;
            //                //cout << EleTag_j - 1 << ",  " << FracTag_j << endl;
            //
            //                _DataType_ Rotate2DTo3D[3][3];
            //                Frac_verts_host[FracTag_j].RoationMatrix(Rotate2DTo3D, 23);
            //
            //                Pos = cuDFNsys::ProductSquare3Vector3<_DataType_>(Rotate2DTo3D, Pos);
            //                Pos = cuDFNsys::MakeVector3(Pos.x + Frac_verts_host[FracTag_j].Center.x,
            //                                            Pos.y + Frac_verts_host[FracTag_j].Center.y,
            //                                            Pos.z + Frac_verts_host[FracTag_j].Center.z);
            //
            //                Position3D[j] = Pos.x;
            //                Position3D[j + numParticles] = Pos.y;
            //                Position3D[j + 2 * numParticles] = Pos.z;
            //                //cout << Pos.x << ", " << Pos.y << ", " << Pos.z << endl;
            //}
            //cout << 2 << endl;
            uint2 dim_pos = make_uint2(3, numParticles);

            h5g.AddDataset(outfilename, "N", "Step_" + cuDFNsys::ToStringWithWidth(i, 10), Position3D.data(),
                           dim_pos);

            //delete[] Position3D;
            //Position3D = NULL;
        };
        //Nproc += 0;
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