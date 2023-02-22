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

#include "ParticleTransport/Transform2DTo3DKernel.cuh"

// ====================================================
// NAME:        Transform2DTo3DKernel
// DESCRIPTION: Transform 2D particle positions to 3D.
// AUTHOR:      Tingchang YIN
// DATE:        22/02/2023
// ====================================================
template <typename T>
__global__ void cuDFNsys::Transform2DTo3DKernel(cuDFNsys::Fracture<T> *Frac_verts_device_ptr,
                                                T *Position3D_dev_ptr,
                                                cuDFNsys::Particle<T> *temp2Dpos_dev_ptr,
                                                uint *ElementFracTag_cuda_devptr,
                                                //uint *EleTag_device_ptr,
                                                uint count)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;

    if (idx > count - 1)
        return;

    uint numParticles = count;
    cuDFNsys::Vector3<T> Pos = cuDFNsys::MakeVector3(temp2Dpos_dev_ptr[idx].Position2D.x,
                                                     temp2Dpos_dev_ptr[idx].Position2D.y,
                                                     (T)0.0);

    //if (idx == 0)
    //printf("temp2Dpos_dev_ptr: %f, %f\n", Pos.x, Pos.y);
    uint EleTag_j = temp2Dpos_dev_ptr[idx].ElementID; //EleTag_device_ptr[idx];

    uint FracTag_j = ElementFracTag_cuda_devptr[EleTag_j - 1];

    //if (idx == 0)
    // printf("%d, FracTag_j: %d\n", idx, FracTag_j);

    T Rotate2DTo3D[3][3];
    Frac_verts_device_ptr[FracTag_j].RoationMatrix(Rotate2DTo3D, 23);

    Pos = cuDFNsys::ProductSquare3Vector3<T>(Rotate2DTo3D, Pos);
    Pos = cuDFNsys::MakeVector3(Pos.x + Frac_verts_device_ptr[FracTag_j].Center.x,
                                Pos.y + Frac_verts_device_ptr[FracTag_j].Center.y,
                                Pos.z + Frac_verts_device_ptr[FracTag_j].Center.z);

    Position3D_dev_ptr[idx] = Pos.x;
    Position3D_dev_ptr[idx + numParticles] = Pos.y;
    Position3D_dev_ptr[idx + 2 * numParticles] = Pos.z;
}; // Transform2DTo3DKernel
template __global__ void cuDFNsys::Transform2DTo3DKernel<double>(cuDFNsys::Fracture<double> *Frac_verts_device_ptr,
                                                                 double *Position3D_dev_ptr,
                                                                 cuDFNsys::Particle<double> *temp2Dpos_dev_ptr,
                                                                 uint *ElementFracTag_cuda_devptr,
                                                                 //uint *EleTag_device_ptr,
                                                                 uint count);
template __global__ void cuDFNsys::Transform2DTo3DKernel<float>(cuDFNsys::Fracture<float> *Frac_verts_device_ptr,
                                                                float *Position3D_dev_ptr,
                                                                cuDFNsys::Particle<float> *temp2Dpos_dev_ptr,
                                                                uint *ElementFracTag_cuda_devptr,
                                                                //uint *EleTag_device_ptr,
                                                                uint count);