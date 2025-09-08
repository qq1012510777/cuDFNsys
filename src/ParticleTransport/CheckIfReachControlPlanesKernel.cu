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

#include "ParticleTransport/CheckIfReachControlPlanesKernel.cuh"

// ====================================================
// NAME:        CheckIfReachControlPlanesKernel
// DESCRIPTION: check if particles reach control planes
// AUTHOR:      Tingchang YIN
// DATE:        11/04/2024
// ====================================================
template <typename T>
__global__ void cuDFNsys::CheckIfReachControlPlanesKernel(
    cuDFNsys::Fracture<T> *Frac_verts_device_ptr,
    cuDFNsys::Particle<T> *temp2Dpos_dev_ptr, uint *ElementFracTag_cuda_devptr,
    //uint *EleTag_device_ptr,
    uint count, uint Dir, T L_percoDir, uint *TimeReachControlPlanes_dev_ptr,
    uint NumControlPlanes, T *ControlPlanes, uint NumParticlesTotal,
    uint StepNo, T *x_ptr, T *y_ptr, T *z_ptr, T magic_number, T InjectionPlaneAtLongitudinalDirection_)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;

    if (idx > count - 1)
        return;

    cuDFNsys::Vector3<T> Pos =
        cuDFNsys::MakeVector3(temp2Dpos_dev_ptr[idx].Position2D.x,
                              temp2Dpos_dev_ptr[idx].Position2D.y, (T)0.0);

    //if (idx == 0)
    //printf("temp2Dpos_dev_ptr: %f, %f\n", Pos.x, Pos.y);
    uint EleTag_j = temp2Dpos_dev_ptr[idx].ElementID; //EleTag_device_ptr[idx];

    uint FracTag_j = ElementFracTag_cuda_devptr[EleTag_j - 1];

    //if (idx == 0)
    // printf("%d, FracTag_j: %d\n", idx, FracTag_j);

    T Rotate2DTo3D[3][3];
    Frac_verts_device_ptr[FracTag_j].RoationMatrix(Rotate2DTo3D, 23);

    Pos = cuDFNsys::ProductSquare3Vector3<T>(Rotate2DTo3D, Pos);
    Pos = cuDFNsys::MakeVector3(
        Pos.x + Frac_verts_device_ptr[FracTag_j].Center.x,
        Pos.y + Frac_verts_device_ptr[FracTag_j].Center.y,
        Pos.z + Frac_verts_device_ptr[FracTag_j].Center.z);

    T *kk = &(Pos.x);
    // kk[Dir] -= (L_percoDir)*temp2Dpos_dev_ptr[idx].FactorPeriodic;

    if (temp2Dpos_dev_ptr[idx].ParticleID > 0 &&
        temp2Dpos_dev_ptr[idx].ParticleID <= NumParticlesTotal)
        x_ptr[idx] = kk[0], y_ptr[idx] = kk[1], z_ptr[idx] = kk[2];
    else
        x_ptr[idx] = magic_number, y_ptr[idx] = magic_number,
        z_ptr[idx] = magic_number;

    for (uint i = 0; i < NumControlPlanes; ++i)
    {
        if (temp2Dpos_dev_ptr[idx].ParticleID > 0 &&
            temp2Dpos_dev_ptr[idx].ParticleID <= NumParticlesTotal &&
            i < NumControlPlanes - 1)
        {
            if (TimeReachControlPlanes_dev_ptr
                    [i * NumParticlesTotal + temp2Dpos_dev_ptr[idx].ParticleID -
                     1] == 0)
            {
                if ((kk[Dir] <= ControlPlanes[i] && ControlPlanes[i] < InjectionPlaneAtLongitudinalDirection_) ||
                    (kk[Dir] >= ControlPlanes[i] && ControlPlanes[i] > InjectionPlaneAtLongitudinalDirection_)
                    )
                {
                    TimeReachControlPlanes_dev_ptr
                        [i * NumParticlesTotal +
                         temp2Dpos_dev_ptr[idx].ParticleID - 1] = StepNo;
                }
                // else if (kk[Dir] <= ControlPlanes[i] && ControlPlanes[i] > InjectionPlaneAtLongitudinalDirection_)
            }
        }
        if (i == NumControlPlanes - 1 &&
            temp2Dpos_dev_ptr[idx].ParticleID < 0 &&
            TimeReachControlPlanes_dev_ptr
                    [i * NumParticlesTotal +
                     abs(temp2Dpos_dev_ptr[idx].ParticleID) - 1] == 0)
        {
            TimeReachControlPlanes_dev_ptr
                [i * NumParticlesTotal +
                 abs(temp2Dpos_dev_ptr[idx].ParticleID) - 1] = StepNo;
        }
    }

}; // check if particles reach control planes
template __global__ void cuDFNsys::CheckIfReachControlPlanesKernel<double>(
    cuDFNsys::Fracture<double> *Frac_verts_device_ptr,
    cuDFNsys::Particle<double> *temp2Dpos_dev_ptr,
    uint *ElementFracTag_cuda_devptr,
    //uint *EleTag_device_ptr,
    uint count, uint Dir, double L_percoDir,
    uint *TimeReachControlPlanes_dev_ptr, uint NumControlPlanes,
    double *ControlPlanes, uint NumParticlesTotal, uint StepNo, double *x_ptr,
    double *y_ptr, double *z_ptr, double magic_number, double InjectionPlaneAtLongitudinalDirection_);
template __global__ void cuDFNsys::CheckIfReachControlPlanesKernel<float>(
    cuDFNsys::Fracture<float> *Frac_verts_device_ptr,
    cuDFNsys::Particle<float> *temp2Dpos_dev_ptr,
    uint *ElementFracTag_cuda_devptr,
    //uint *EleTag_device_ptr,
    uint count, uint Dir, float L_percoDir,
    uint *TimeReachControlPlanes_dev_ptr, uint NumControlPlanes,
    float *ControlPlanes, uint NumParticlesTotal, uint StepNo, float *x_ptr,
    float *y_ptr, float *z_ptr, float magic_number, float InjectionPlaneAtLongitudinalDirection_);