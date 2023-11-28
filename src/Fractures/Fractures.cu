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

#include "Fractures/Fractures.cuh"

// ====================================================
// NAME:        Fractures
// DESCRIPTION: Fractures in a DFN
// AUTHOR:      Tingchang YIN
// DATE:        04/04/2022
// ====================================================
template <typename T>
__global__ void cuDFNsys::Fractures(
    cuDFNsys::Fracture<T> *verts, unsigned long seed, int count,
    cuDFNsys::Vector1<T> model_L,
    uint
        ModeSizeDistri, // 0 = power law; 1 = lognormal; 2 = uniform; 3 = monosize
    cuDFNsys::Vector4<T> ParaSizeDistri, // when mode = 1, ;
    cuDFNsys::Vector1<T> kappa,
    cuDFNsys::Vector1<T> conductivity_powerlaw_exponent, T Gamma_constant,
    double3 DimensionRatio, cuDFNsys::Vector3<T> MeanOrientation)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;

    if (i > count - 1)
        return;

    curandState state;

    curand_init(seed, i, 0, &state);

    cuDFNsys::Vector1<T> R_ = 0;

    // if (alpha == 0 && abs(minR - maxR) < 1e-7)
    //     R_ = minR;
    // else if (alpha == 0 && abs(minR - maxR) > 1e-7)
    //     R_ = cuDFNsys::RandomUniform(minR, maxR, curand_uniform(&state));
    // else
    //     R_ = cuDFNsys::RandomPowerlaw(minR, maxR, alpha, curand_uniform(&state));

    if (ModeSizeDistri == 0)
        R_ = cuDFNsys::RandomPowerlaw((T)ParaSizeDistri.y, (T)ParaSizeDistri.z,
                                      (T)ParaSizeDistri.x,
                                      (T)curand_uniform(&state));
    else if (ModeSizeDistri == 1)
        R_ = cuDFNsys::RandomLognormal((T)ParaSizeDistri.x, (T)ParaSizeDistri.y,
                                       (T)ParaSizeDistri.z, (T)ParaSizeDistri.w,
                                       (T)curand_uniform(&state));
    else if (ModeSizeDistri == 2)
        R_ = cuDFNsys::RandomUniform((T)ParaSizeDistri.x, (T)ParaSizeDistri.y,
                                     (T)curand_uniform(&state));
    else if (ModeSizeDistri == 3)
        R_ = ParaSizeDistri.x;

    verts[i].Radius = R_;
    //printf("%f\n", verts[i].Radius);

    verts[i].NumVertsTruncated = 4;
    //printf("conductivity_powerlaw_exponent: %f, kappa: %f\n", conductivity_powerlaw_exponent, kappa);

    if (conductivity_powerlaw_exponent == 0)
        verts[i].Conductivity = (T)(pow(1.0e-3, 3.0) / 12.0);
    else
    {
        verts[i].Conductivity = (pow(Gamma_constant, 3.0) / 12.0) *
                                pow(R_, 3.0 * conductivity_powerlaw_exponent);
    }

    if (DimensionRatio.x == 1 && DimensionRatio.y == 1 && DimensionRatio.z == 1)
    {
        verts[i].Center.x =
            cuDFNsys::RandomUniform((T)(-model_L * 0.5), (T)(model_L * 0.5),
                                    (T)(curand_uniform(&state)));
        verts[i].Center.y =
            cuDFNsys::RandomUniform((T)(-model_L * 0.5), (T)(model_L * 0.5),
                                    (T)(curand_uniform(&state)));
        verts[i].Center.z =
            cuDFNsys::RandomUniform((T)(-model_L * 0.5), (T)(model_L * 0.5),
                                    (T)(curand_uniform(&state)));
    }
    else
    {
        verts[i].Center.x = cuDFNsys::RandomUniform(
            (T)(-DimensionRatio.x * model_L * 0.5),
            (T)(DimensionRatio.x * model_L * 0.5), (T)(curand_uniform(&state)));
        verts[i].Center.y = cuDFNsys::RandomUniform(
            (T)(-DimensionRatio.y * model_L * 0.5),
            (T)(DimensionRatio.y * model_L * 0.5), (T)(curand_uniform(&state)));
        verts[i].Center.z = cuDFNsys::RandomUniform(
            (T)(-DimensionRatio.z * model_L * 0.5),
            (T)(DimensionRatio.z * model_L * 0.5), (T)(curand_uniform(&state)));
    }

    verts[i].NormalVec = cuDFNsys::MakeVector3(
        (T)cuDFNsys::RandomUniform(
            (cuDFNsys::Vector1<T>)-1.0, (cuDFNsys::Vector1<T>)1.0,
            (cuDFNsys::Vector1<T>)curand_uniform(&state)),
        (T)cuDFNsys::RandomUniform(
            (cuDFNsys::Vector1<T>)-1.0, (cuDFNsys::Vector1<T>)1.0,
            (cuDFNsys::Vector1<T>)curand_uniform(&state)),
        (T)0);
    cuDFNsys::Vector1<T> R_xy =
        sqrt(verts[i].NormalVec.x * verts[i].NormalVec.x +
             verts[i].NormalVec.y * verts[i].NormalVec.y);
    verts[i].NormalVec.z =
        R_xy / tan(cuDFNsys::RandomFisher((T)curand_uniform(&state), (T)kappa));
    cuDFNsys::Vector1<T> norm_f =
        sqrt(verts[i].NormalVec.x * verts[i].NormalVec.x +
             verts[i].NormalVec.y * verts[i].NormalVec.y +
             verts[i].NormalVec.z * verts[i].NormalVec.z);
    verts[i].NormalVec.x /= norm_f;
    verts[i].NormalVec.y /= norm_f;
    verts[i].NormalVec.z /= norm_f;

    //---------------------- if the mean orientation is not (0, 0, 1)
    if (!(MeanOrientation.x == (T)0. && MeanOrientation.y == (T)0. &&
          MeanOrientation.z == (T)1.))
    {
        //rotate the orientation of this fracture
        T angle_fs = acos(MeanOrientation.z);
        cuDFNsys::Vector3<T> Ori_Rotate = cuDFNsys::CrossProductVector3<T>(
            cuDFNsys::MakeVector3((T)0., (T)0., (T)1.), MeanOrientation);

        T norm_ori =
            pow(Ori_Rotate.x * Ori_Rotate.x + Ori_Rotate.y * Ori_Rotate.y +
                    Ori_Rotate.z * Ori_Rotate.z,
                0.5);
        Ori_Rotate.x /= norm_ori, Ori_Rotate.y /= norm_ori,
            Ori_Rotate.z /= norm_ori;

        cuDFNsys::Quaternion<T> qua;
        qua = qua.DescribeRotation(Ori_Rotate, angle_fs);

        verts[i].NormalVec = qua.Rotate(verts[i].NormalVec);

        norm_f = sqrt(verts[i].NormalVec.x * verts[i].NormalVec.x +
                      verts[i].NormalVec.y * verts[i].NormalVec.y +
                      verts[i].NormalVec.z * verts[i].NormalVec.z);
        verts[i].NormalVec.x /= norm_f;
        verts[i].NormalVec.y /= norm_f;
        verts[i].NormalVec.z /= norm_f;

        if (verts[i].NormalVec.z < 0)
            verts[i].NormalVec.x *= -1, verts[i].NormalVec.y *= -1,
                verts[i].NormalVec.z *= -1;
    }

    cuDFNsys::Vector1<T> *normal_fff = &verts[i].NormalVec.x;
    cuDFNsys::Vector1<T> *verts_3D_ptr = &verts[i].Verts3D[0].x;
    for (int j = 0; j < 3; ++j)
    {
        if (abs(normal_fff[j]) > 1e-3)
        {
            verts_3D_ptr[(j + 1) % 3] = cuDFNsys::RandomUniform(
                (cuDFNsys::Vector1<T>)-1.0, (cuDFNsys::Vector1<T>)1.0,
                (cuDFNsys::Vector1<T>)curand_uniform(&state));
            verts_3D_ptr[(j + 2) % 3] = cuDFNsys::RandomUniform(
                (cuDFNsys::Vector1<T>)-1.0, (cuDFNsys::Vector1<T>)1.0,
                (cuDFNsys::Vector1<T>)curand_uniform(&state));
            verts_3D_ptr[j] =
                -1.0 *
                (verts_3D_ptr[(j + 1) % 3] * normal_fff[(j + 1) % 3] +
                 verts_3D_ptr[(j + 2) % 3] * normal_fff[(j + 2) % 3]) /
                normal_fff[j];
            break;
        }
    }

    cuDFNsys::Vector1<T> norm_vert1 =
        sqrt(verts[i].Verts3D[0].x * verts[i].Verts3D[0].x +
             verts[i].Verts3D[0].y * verts[i].Verts3D[0].y +
             verts[i].Verts3D[0].z * verts[i].Verts3D[0].z);
    norm_vert1 = R_ / norm_vert1;
    verts[i].Verts3D[0].x *= norm_vert1;
    verts[i].Verts3D[0].y *= norm_vert1;
    verts[i].Verts3D[0].z *= norm_vert1;
    verts[i].Verts3D[2].x = -1.0 * verts[i].Verts3D[0].x;
    verts[i].Verts3D[2].y = -1.0 * verts[i].Verts3D[0].y;
    verts[i].Verts3D[2].z = -1.0 * verts[i].Verts3D[0].z;

    verts[i].Verts3D[1] = cuDFNsys::CrossProductVector3<T>(verts[i].NormalVec,
                                                           verts[i].Verts3D[0]);
    norm_vert1 = sqrt(verts[i].Verts3D[1].x * verts[i].Verts3D[1].x +
                      verts[i].Verts3D[1].y * verts[i].Verts3D[1].y +
                      verts[i].Verts3D[1].z * verts[i].Verts3D[1].z);
    norm_vert1 = R_ / norm_vert1;
    verts[i].Verts3D[1].x *= norm_vert1;
    verts[i].Verts3D[1].y *= norm_vert1;
    verts[i].Verts3D[1].z *= norm_vert1;
    verts[i].Verts3D[3].x = -1.0 * verts[i].Verts3D[1].x;
    verts[i].Verts3D[3].y = -1.0 * verts[i].Verts3D[1].y;
    verts[i].Verts3D[3].z = -1.0 * verts[i].Verts3D[1].z;
    //-----------------------------------------
    for (int j = 0; j < 4; ++j)
    {
        verts[i].Verts3D[j].x += verts[i].Center.x;
        verts[i].Verts3D[j].y += verts[i].Center.y;
        verts[i].Verts3D[j].z += verts[i].Center.z;

        verts[i].Verts3DTruncated[j].x = verts[i].Verts3D[j].x;
        verts[i].Verts3DTruncated[j].y = verts[i].Verts3D[j].y;
        verts[i].Verts3DTruncated[j].z = verts[i].Verts3D[j].z;
    };

    bool gh = cuDFNsys::TruncateFracture<T>(&verts[i],
                                            DimensionRatio.x * model_L, 0, -1);
    verts[i].ConnectModelSurf[0] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], DimensionRatio.x * model_L, 0,
                                       1);
    verts[i].ConnectModelSurf[1] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], DimensionRatio.y * model_L, 1,
                                       -1);
    verts[i].ConnectModelSurf[2] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], DimensionRatio.y * model_L, 1,
                                       1);
    verts[i].ConnectModelSurf[3] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], DimensionRatio.z * model_L, 2,
                                       -1);
    verts[i].ConnectModelSurf[4] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], DimensionRatio.z * model_L, 2,
                                       1);
    verts[i].ConnectModelSurf[5] = gh;
}; // Fractures
template __global__ void cuDFNsys::Fractures<double>(
    cuDFNsys::Fracture<double> *verts, unsigned long seed, int count,
    cuDFNsys::Vector1<double> model_L,
    uint
        ModeSizeDistri, // 1 = power law; 2 = lognormal; 3 = uniform; 4 = monosize
    cuDFNsys::Vector4<double> ParaSizeDistri, // when mode = 1, ;
    cuDFNsys::Vector1<double> kappa,
    cuDFNsys::Vector1<double> conductivity_powerlaw_exponent,
    double Gamma_constant, double3 DimensionRatio,
    cuDFNsys::Vector3<double> MeanOrientation);
template __global__ void cuDFNsys::Fractures<float>(
    cuDFNsys::Fracture<float> *verts, unsigned long seed, int count,
    cuDFNsys::Vector1<float> model_L,
    uint
        ModeSizeDistri, // 1 = power law; 2 = lognormal; 3 = uniform; 4 = monosize
    cuDFNsys::Vector4<float> ParaSizeDistri, // when mode = 1, ;
    cuDFNsys::Vector1<float> kappa,
    cuDFNsys::Vector1<float> conductivity_powerlaw_exponent,
    float Gamma_constant, double3 DimensionRatio,
    cuDFNsys::Vector3<float> MeanOrientation);

// ====================================================
// NAME:        FracturesDeterministic
// DESCRIPTION: Fractures in a DFN
// AUTHOR:      Tingchang YIN
// DATE:        03/11/2023
// ====================================================
template <typename T>
__global__ void cuDFNsys::FracturesDeterministic<T>(
    cuDFNsys::Fracture<T> *verts, unsigned long seed, T *Data_f, int count,
    T model_L, int PercoDir, double3 DomainDimensionRatio)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;

    if (i > count - 1)
        return;

    int DataNumForOneFracture = 9;

    verts[i].NormalVec.x = Data_f[i * DataNumForOneFracture],
    verts[i].NormalVec.y = Data_f[i * DataNumForOneFracture + 1],
    verts[i].NormalVec.z = Data_f[i * DataNumForOneFracture + 2],
    verts[i].Center.x = Data_f[i * DataNumForOneFracture + 3],
    verts[i].Center.y = Data_f[i * DataNumForOneFracture + 4],
    verts[i].Center.z = Data_f[i * DataNumForOneFracture + 5],
    verts[i].Radius = Data_f[i * DataNumForOneFracture + 6];

    T conductivity_powerlaw_exponent = Data_f[i * DataNumForOneFracture + 7],
      Gamma_constant = Data_f[i * DataNumForOneFracture + 8];

    cuDFNsys::Vector1<T> norm_f =
        sqrt(verts[i].NormalVec.x * verts[i].NormalVec.x +
             verts[i].NormalVec.y * verts[i].NormalVec.y +
             verts[i].NormalVec.z * verts[i].NormalVec.z);
    verts[i].NormalVec.x /= norm_f;
    verts[i].NormalVec.y /= norm_f;
    verts[i].NormalVec.z /= norm_f;

    if (verts[i].NormalVec.z < 0)
        verts[i].NormalVec.x *= -1, verts[i].NormalVec.y *= -1,
            verts[i].NormalVec.z *= -1;

    //---------------------

    verts[i].NumVertsTruncated = 4;
    if (conductivity_powerlaw_exponent == 0)
        verts[i].Conductivity = (T)(pow(1.0e-3, 3.0) / 12.0);
    else
    {
        verts[i].Conductivity =
            (pow(Gamma_constant, 3.0) / 12.0) *
            pow(verts[i].Radius, 3.0 * conductivity_powerlaw_exponent);
    }
    //----------------
    curandState state;

    curand_init(seed, i, 0, &state);

    // generate vertices
    cuDFNsys::Vector1<T> *normal_fff = &verts[i].NormalVec.x;
    cuDFNsys::Vector1<T> *verts_3D_ptr = &verts[i].Verts3D[0].x;
    for (int j = 0; j < 3; ++j)
    {
        if (abs(normal_fff[j]) > 1e-3)
        {
            verts_3D_ptr[(j + 1) % 3] = cuDFNsys::RandomUniform(
                (cuDFNsys::Vector1<T>)-1.0, (cuDFNsys::Vector1<T>)1.0,
                (cuDFNsys::Vector1<T>)curand_uniform(&state));
            verts_3D_ptr[(j + 2) % 3] = cuDFNsys::RandomUniform(
                (cuDFNsys::Vector1<T>)-1.0, (cuDFNsys::Vector1<T>)1.0,
                (cuDFNsys::Vector1<T>)curand_uniform(&state));
            verts_3D_ptr[j] =
                -1.0 *
                (verts_3D_ptr[(j + 1) % 3] * normal_fff[(j + 1) % 3] +
                 verts_3D_ptr[(j + 2) % 3] * normal_fff[(j + 2) % 3]) /
                normal_fff[j];
            break;
        }
    }

    cuDFNsys::Vector1<T> norm_vert1 =
        sqrt(verts[i].Verts3D[0].x * verts[i].Verts3D[0].x +
             verts[i].Verts3D[0].y * verts[i].Verts3D[0].y +
             verts[i].Verts3D[0].z * verts[i].Verts3D[0].z);
    norm_vert1 = verts[i].Radius / norm_vert1;
    verts[i].Verts3D[0].x *= norm_vert1;
    verts[i].Verts3D[0].y *= norm_vert1;
    verts[i].Verts3D[0].z *= norm_vert1;
    verts[i].Verts3D[2].x = -1.0 * verts[i].Verts3D[0].x;
    verts[i].Verts3D[2].y = -1.0 * verts[i].Verts3D[0].y;
    verts[i].Verts3D[2].z = -1.0 * verts[i].Verts3D[0].z;

    verts[i].Verts3D[1] = cuDFNsys::CrossProductVector3<T>(verts[i].NormalVec,
                                                           verts[i].Verts3D[0]);
    norm_vert1 = sqrt(verts[i].Verts3D[1].x * verts[i].Verts3D[1].x +
                      verts[i].Verts3D[1].y * verts[i].Verts3D[1].y +
                      verts[i].Verts3D[1].z * verts[i].Verts3D[1].z);
    norm_vert1 = verts[i].Radius / norm_vert1;
    verts[i].Verts3D[1].x *= norm_vert1;
    verts[i].Verts3D[1].y *= norm_vert1;
    verts[i].Verts3D[1].z *= norm_vert1;
    verts[i].Verts3D[3].x = -1.0 * verts[i].Verts3D[1].x;
    verts[i].Verts3D[3].y = -1.0 * verts[i].Verts3D[1].y;
    verts[i].Verts3D[3].z = -1.0 * verts[i].Verts3D[1].z;
    //---------------
    //-----------------------------------------
    for (int j = 0; j < 4; ++j)
    {
        verts[i].Verts3D[j].x += verts[i].Center.x;
        verts[i].Verts3D[j].y += verts[i].Center.y;
        verts[i].Verts3D[j].z += verts[i].Center.z;

        verts[i].Verts3DTruncated[j].x = verts[i].Verts3D[j].x;
        verts[i].Verts3DTruncated[j].y = verts[i].Verts3D[j].y;
        verts[i].Verts3DTruncated[j].z = verts[i].Verts3D[j].z;
    };
    //-----------------------
    bool gh = cuDFNsys::TruncateFracture<T>(
        &verts[i], DomainDimensionRatio.x * model_L, 0, -1);
    verts[i].ConnectModelSurf[0] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i],
                                       DomainDimensionRatio.x * model_L, 0, 1);
    verts[i].ConnectModelSurf[1] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i],
                                       DomainDimensionRatio.y * model_L, 1, -1);
    verts[i].ConnectModelSurf[2] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i],
                                       DomainDimensionRatio.y * model_L, 1, 1);
    verts[i].ConnectModelSurf[3] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i],
                                       DomainDimensionRatio.z * model_L, 2, -1);
    verts[i].ConnectModelSurf[4] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i],
                                       DomainDimensionRatio.z * model_L, 2, 1);
    verts[i].ConnectModelSurf[5] = gh;
};
template __global__ void cuDFNsys::FracturesDeterministic<double>(
    cuDFNsys::Fracture<double> *verts, unsigned long seed, double *Data_f,
    int count, double model_L, int PercoDir, double3 DomainDimensionRatio);
template __global__ void cuDFNsys::FracturesDeterministic<float>(
    cuDFNsys::Fracture<float> *verts, unsigned long seed, float *Data_f,
    int count, float model_L, int PercoDir, double3 DomainDimensionRatio);

// ====================================================
// NAME:        FracturesForSpatialPeriodicity
// DESCRIPTION: Fractures in a DFN
// AUTHOR:      Tingchang YIN
// DATE:        06/11/2023
// ====================================================
template <typename T>
__global__ void cuDFNsys::FracturesForSpatialPeriodicity<T>(
    cuDFNsys::Fracture<T> *verts, int count, cuDFNsys::Vector1<T> model_L,
    int PercoDir, double3 DomainDimensionRatio)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;

    if (i > count - 1)
        return;

    T *Center_fff = &verts[i].Center.x;
    double *Dimensionratio = &DomainDimensionRatio.x;

    T Increament_ff = model_L * Dimensionratio[PercoDir] / 2;

    Center_fff[PercoDir] += Increament_ff;
    verts[i].NumVertsTruncated = 4;

    for (int j = 0; j < 4; ++j)
    {
        T *Vertic_ff = &(verts[i].Verts3D[j].x);
        Vertic_ff[PercoDir] += Increament_ff;

        verts[i].Verts3DTruncated[j].x = verts[i].Verts3D[j].x,
        verts[i].Verts3DTruncated[j].y = verts[i].Verts3D[j].y,
        verts[i].Verts3DTruncated[j].z = verts[i].Verts3D[j].z;
    }

    if (i >= count / 2)
    {
        Center_fff[PercoDir] = -Center_fff[PercoDir];

        for (int j = 0; j < 4; ++j)
        {
            T *Vertic_ff = &(verts[i].Verts3D[j].x);
            Vertic_ff[PercoDir] = -Vertic_ff[PercoDir];

            verts[i].Verts3DTruncated[j].x = verts[i].Verts3D[j].x,
            verts[i].Verts3DTruncated[j].y = verts[i].Verts3D[j].y,
            verts[i].Verts3DTruncated[j].z = verts[i].Verts3D[j].z;
        }

        cuDFNsys::Vector3<T> D1 = cuDFNsys::MakeVector3(
                                 verts[i].Verts3D[0].x - verts[i].Center.x,
                                 verts[i].Verts3D[0].y - verts[i].Center.y,
                                 verts[i].Verts3D[0].z - verts[i].Center.z),
                             D2 = cuDFNsys::MakeVector3(
                                 verts[i].Verts3D[1].x - verts[i].Center.x,
                                 verts[i].Verts3D[1].y - verts[i].Center.y,
                                 verts[i].Verts3D[1].z - verts[i].Center.z);
        verts[i].NormalVec = cuDFNsys::CrossProductVector3<T>(D1, D2);

        cuDFNsys::Vector1<T> norm_f =
            sqrt(verts[i].NormalVec.x * verts[i].NormalVec.x +
                 verts[i].NormalVec.y * verts[i].NormalVec.y +
                 verts[i].NormalVec.z * verts[i].NormalVec.z);

        verts[i].NormalVec.x /= norm_f;
        verts[i].NormalVec.y /= norm_f;
        verts[i].NormalVec.z /= norm_f;

        if (verts[i].NormalVec.z < 0)
            verts[i].NormalVec.x = -verts[i].NormalVec.x,
            verts[i].NormalVec.y = -verts[i].NormalVec.y,
            verts[i].NormalVec.z = -verts[i].NormalVec.z;
    }

    verts[i].ConnectModelSurf[0] = false, verts[i].ConnectModelSurf[1] = false,
    verts[i].ConnectModelSurf[2] = false, verts[i].ConnectModelSurf[3] = false,
    verts[i].ConnectModelSurf[4] = false, verts[i].ConnectModelSurf[5] = false;

    cuDFNsys::Vector3<T> DomainElongate =
        cuDFNsys::MakeVector3((T)1., (T)1., (T)1.);
    T *DomainElongate_ff = &(DomainElongate.x);
    DomainElongate_ff[PercoDir] = 2;

    bool gh = cuDFNsys::TruncateFracture<T>(
        &verts[i], DomainDimensionRatio.x * model_L * DomainElongate_ff[0], 0,
        -1);
    verts[i].ConnectModelSurf[0] = gh;

    gh = cuDFNsys::TruncateFracture<T>(
        &verts[i], DomainDimensionRatio.x * model_L * DomainElongate_ff[0], 0,
        1);
    verts[i].ConnectModelSurf[1] = gh;

    gh = cuDFNsys::TruncateFracture<T>(
        &verts[i], DomainDimensionRatio.y * model_L * DomainElongate_ff[1], 1,
        -1);
    verts[i].ConnectModelSurf[2] = gh;

    gh = cuDFNsys::TruncateFracture<T>(
        &verts[i], DomainDimensionRatio.y * model_L * DomainElongate_ff[1], 1,
        1);
    verts[i].ConnectModelSurf[3] = gh;

    gh = cuDFNsys::TruncateFracture<T>(
        &verts[i], DomainDimensionRatio.z * model_L * DomainElongate_ff[2], 2,
        -1);
    verts[i].ConnectModelSurf[4] = gh;

    gh = cuDFNsys::TruncateFracture<T>(
        &verts[i], DomainDimensionRatio.z * model_L * DomainElongate_ff[2], 2,
        1);
    verts[i].ConnectModelSurf[5] = gh;
};
template __global__ void cuDFNsys::FracturesForSpatialPeriodicity<double>(
    cuDFNsys::Fracture<double> *verts, int count,
    cuDFNsys::Vector1<double> model_L, int PercoDir,
    double3 DomainDimensionRatio);
template __global__ void cuDFNsys::FracturesForSpatialPeriodicity<float>(
    cuDFNsys::Fracture<float> *verts, int count,
    cuDFNsys::Vector1<float> model_L, int PercoDir,
    double3 DomainDimensionRatio);

// ====================================================
// NAME:        Fractures
// DESCRIPTION: Fractures in a DFN
// AUTHOR:      Tingchang YIN
// DATE:        04/04/2022
// ====================================================
template <typename T>
cuDFNsys::FracturesCPU<T>::FracturesCPU(
    thrust::host_vector<cuDFNsys::Fracture<T>> &verts, unsigned long seed,
    int count, T model_L,
    uint
        ModeSizeDistri, // 0 = power law; 1 = lognormal; 2 = uniform; 3 = monosize
    cuDFNsys::Vector4<T> ParaSizeDistri, // when mode = 1, ;
    T kappa, T conductivity_powerlaw_exponent, uint Nproc)
{
    srand(seed);
    throw cuDFNsys::ExceptionsPause(" `FracturesCPU` is Not recommended!\n");
#pragma omp parallel for schedule(static) num_threads(Nproc)
    for (int i = 0; i < count; ++i)
    {
        T R_ = 0;

        if (ModeSizeDistri == 0)
            R_ = cuDFNsys::RandomPowerlaw(
                (T)ParaSizeDistri.y, (T)ParaSizeDistri.z, (T)ParaSizeDistri.x,
                (T)((T)rand() / (T)RAND_MAX));
        else if (ModeSizeDistri == 1)
            R_ = cuDFNsys::RandomLognormal(
                (T)ParaSizeDistri.x, (T)ParaSizeDistri.y, (T)ParaSizeDistri.z,
                (T)ParaSizeDistri.w, (T)((T)rand() / (T)RAND_MAX));
        else if (ModeSizeDistri == 2)
            R_ = cuDFNsys::RandomUniform((T)ParaSizeDistri.x,
                                         (T)ParaSizeDistri.y,
                                         (T)((T)rand() / (T)RAND_MAX));
        else if (ModeSizeDistri == 3)
            R_ = ParaSizeDistri.x;

        verts[i].Radius = R_;
        //printf("%f\n", verts[i].Radius);

        verts[i].NumVertsTruncated = 4;
        //printf("conductivity_powerlaw_exponent: %f, kappa: %f\n", conductivity_powerlaw_exponent, kappa);

        if (conductivity_powerlaw_exponent == 0)
            verts[i].Conductivity = (T)(pow(1.0e-3, 3.0) / 12.0);
        else
            verts[i].Conductivity =
                (1.0e-11) * pow(R_, 3.0 * conductivity_powerlaw_exponent);

        verts[i].Center.x = cuDFNsys::RandomUniform(
            (T)(-model_L * 0.5), (T)(model_L * 0.5), ((T)rand() / (T)RAND_MAX));
        verts[i].Center.y = cuDFNsys::RandomUniform(
            (T)(-model_L * 0.5), (T)(model_L * 0.5), ((T)rand() / (T)RAND_MAX));
        verts[i].Center.z = cuDFNsys::RandomUniform(
            (T)(-model_L * 0.5), (T)(model_L * 0.5), ((T)rand() / (T)RAND_MAX));

        verts[i].NormalVec = cuDFNsys::MakeVector3(
            (T)cuDFNsys::RandomUniform((cuDFNsys::Vector1<T>)-1.0,
                                       (cuDFNsys::Vector1<T>)1.0,
                                       ((T)rand() / (T)RAND_MAX)),
            (T)cuDFNsys::RandomUniform((cuDFNsys::Vector1<T>)-1.0,
                                       (cuDFNsys::Vector1<T>)1.0,
                                       ((T)rand() / (T)RAND_MAX)),
            (T)0);
        cuDFNsys::Vector1<T> R_xy =
            sqrt(verts[i].NormalVec.x * verts[i].NormalVec.x +
                 verts[i].NormalVec.y * verts[i].NormalVec.y);
        verts[i].NormalVec.z =
            R_xy /
            tan(cuDFNsys::RandomFisher(((T)rand() / (T)RAND_MAX), (T)kappa));
        cuDFNsys::Vector1<T> norm_f =
            sqrt(verts[i].NormalVec.x * verts[i].NormalVec.x +
                 verts[i].NormalVec.y * verts[i].NormalVec.y +
                 verts[i].NormalVec.z * verts[i].NormalVec.z);
        verts[i].NormalVec.x /= norm_f;
        verts[i].NormalVec.y /= norm_f;
        verts[i].NormalVec.z /= norm_f;

        cuDFNsys::Vector1<T> *normal_fff = &verts[i].NormalVec.x;
        cuDFNsys::Vector1<T> *verts_3D_ptr = &verts[i].Verts3D[0].x;
        for (int j = 0; j < 3; ++j)
        {
            if (abs(normal_fff[j]) > 1e-3)
            {
                verts_3D_ptr[(j + 1) % 3] = cuDFNsys::RandomUniform(
                    (cuDFNsys::Vector1<T>)-1.0, (cuDFNsys::Vector1<T>)1.0,
                    ((T)rand() / (T)RAND_MAX));
                verts_3D_ptr[(j + 2) % 3] = cuDFNsys::RandomUniform(
                    (cuDFNsys::Vector1<T>)-1.0, (cuDFNsys::Vector1<T>)1.0,
                    ((T)rand() / (T)RAND_MAX));
                verts_3D_ptr[j] =
                    -1.0 *
                    (verts_3D_ptr[(j + 1) % 3] * normal_fff[(j + 1) % 3] +
                     verts_3D_ptr[(j + 2) % 3] * normal_fff[(j + 2) % 3]) /
                    normal_fff[j];
                break;
            }
        }

        cuDFNsys::Vector1<T> norm_vert1 =
            sqrt(verts[i].Verts3D[0].x * verts[i].Verts3D[0].x +
                 verts[i].Verts3D[0].y * verts[i].Verts3D[0].y +
                 verts[i].Verts3D[0].z * verts[i].Verts3D[0].z);
        norm_vert1 = R_ / norm_vert1;
        verts[i].Verts3D[0].x *= norm_vert1;
        verts[i].Verts3D[0].y *= norm_vert1;
        verts[i].Verts3D[0].z *= norm_vert1;
        verts[i].Verts3D[2].x = -1.0 * verts[i].Verts3D[0].x;
        verts[i].Verts3D[2].y = -1.0 * verts[i].Verts3D[0].y;
        verts[i].Verts3D[2].z = -1.0 * verts[i].Verts3D[0].z;

        verts[i].Verts3D[1] = cuDFNsys::CrossProductVector3<T>(
            verts[i].NormalVec, verts[i].Verts3D[0]);
        norm_vert1 = sqrt(verts[i].Verts3D[1].x * verts[i].Verts3D[1].x +
                          verts[i].Verts3D[1].y * verts[i].Verts3D[1].y +
                          verts[i].Verts3D[1].z * verts[i].Verts3D[1].z);
        norm_vert1 = R_ / norm_vert1;
        verts[i].Verts3D[1].x *= norm_vert1;
        verts[i].Verts3D[1].y *= norm_vert1;
        verts[i].Verts3D[1].z *= norm_vert1;
        verts[i].Verts3D[3].x = -1.0 * verts[i].Verts3D[1].x;
        verts[i].Verts3D[3].y = -1.0 * verts[i].Verts3D[1].y;
        verts[i].Verts3D[3].z = -1.0 * verts[i].Verts3D[1].z;
        //-----------------------------------------
        for (int j = 0; j < 4; ++j)
        {
            verts[i].Verts3D[j].x += verts[i].Center.x;
            verts[i].Verts3D[j].y += verts[i].Center.y;
            verts[i].Verts3D[j].z += verts[i].Center.z;

            verts[i].Verts3DTruncated[j].x = verts[i].Verts3D[j].x;
            verts[i].Verts3DTruncated[j].y = verts[i].Verts3D[j].y;
            verts[i].Verts3DTruncated[j].z = verts[i].Verts3D[j].z;
        };

        bool gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 0, -1);
        verts[i].ConnectModelSurf[0] = gh;

        gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 0, 1);
        verts[i].ConnectModelSurf[1] = gh;

        gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 1, -1);
        verts[i].ConnectModelSurf[2] = gh;

        gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 1, 1);
        verts[i].ConnectModelSurf[3] = gh;

        gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 2, -1);
        verts[i].ConnectModelSurf[4] = gh;

        gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 2, 1);
        verts[i].ConnectModelSurf[5] = gh;
    }
}; // FracturesCPU
template cuDFNsys::FracturesCPU<double>::FracturesCPU(
    thrust::host_vector<cuDFNsys::Fracture<double>> &verts, unsigned long seed,
    int count, double model_L,
    uint
        ModeSizeDistri, // 0 = power law; 1 = lognormal; 2 = uniform; 3 = monosize
    cuDFNsys::Vector4<double> ParaSizeDistri, // when mode = 1, ;
    double kappa, double conductivity_powerlaw_exponent, uint Nproc);
template cuDFNsys::FracturesCPU<float>::FracturesCPU(
    thrust::host_vector<cuDFNsys::Fracture<float>> &verts, unsigned long seed,
    int count, float model_L,
    uint
        ModeSizeDistri, // 0 = power law; 1 = lognormal; 2 = uniform; 3 = monosize
    cuDFNsys::Vector4<float> ParaSizeDistri, // when mode = 1, ;
    float kappa, float conductivity_powerlaw_exponent, uint Nproc);

// ====================================================
// NAME:        FracturesCrossedVertical
// DESCRIPTION: two crossed vertical fractures
// AUTHOR:      Tingchang YIN
// DATE:        20/04/2022
// ====================================================
template <typename T>
__global__ void cuDFNsys::FracturesCrossedVertical(cuDFNsys::Fracture<T> *verts,
                                                   unsigned long seed,
                                                   int count, T model_L)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;

    if (i >= 2)
        return;

    curandState state;

    curand_init(seed, i, 0, &state);

    T R_ = sqrt(model_L * model_L / 2.0f);

    verts[i].Radius = R_;
    //printf("%f\n", verts[i].Radius);

    verts[i].NumVertsTruncated = 4;
    //printf("conductivity_powerlaw_exponent: %f, kappa: %f\n", conductivity_powerlaw_exponent, kappa);

    verts[i].Conductivity = (T)(pow(1.0e-3, 3.0) / 12.0);

    verts[i].Center.x = 0;
    verts[i].Center.y = 0;
    verts[i].Center.z = 0;

    if (i == 0)
        verts[i].NormalVec = cuDFNsys::MakeVector3((T)0.0, (T)1.0, (T)0.0);
    else
        verts[i].NormalVec = cuDFNsys::MakeVector3((T)1.0, (T)0.0, (T)0.0);

    if (i == 0)
        verts[i].Verts3D[0] = cuDFNsys::MakeVector3((T)(-0.5 * model_L), (T)0.0,
                                                    (T)(0.5 * model_L));
    else
        verts[i].Verts3D[0] = cuDFNsys::MakeVector3((T)0.0, (T)(-0.5 * model_L),
                                                    (T)(0.5 * model_L));

    T norm_vert1 = sqrt(verts[i].Verts3D[0].x * verts[i].Verts3D[0].x +
                        verts[i].Verts3D[0].y * verts[i].Verts3D[0].y +
                        verts[i].Verts3D[0].z * verts[i].Verts3D[0].z);
    norm_vert1 = R_ / norm_vert1;
    verts[i].Verts3D[0].x *= norm_vert1;
    verts[i].Verts3D[0].y *= norm_vert1;
    verts[i].Verts3D[0].z *= norm_vert1;
    verts[i].Verts3D[2].x = -1.0 * verts[i].Verts3D[0].x;
    verts[i].Verts3D[2].y = -1.0 * verts[i].Verts3D[0].y;
    verts[i].Verts3D[2].z = -1.0 * verts[i].Verts3D[0].z;

    verts[i].Verts3D[1] = cuDFNsys::CrossProductVector3<T>(verts[i].NormalVec,
                                                           verts[i].Verts3D[0]);
    norm_vert1 = sqrt(verts[i].Verts3D[1].x * verts[i].Verts3D[1].x +
                      verts[i].Verts3D[1].y * verts[i].Verts3D[1].y +
                      verts[i].Verts3D[1].z * verts[i].Verts3D[1].z);
    norm_vert1 = R_ / norm_vert1;
    verts[i].Verts3D[1].x *= norm_vert1;
    verts[i].Verts3D[1].y *= norm_vert1;
    verts[i].Verts3D[1].z *= norm_vert1;
    verts[i].Verts3D[3].x = -1.0 * verts[i].Verts3D[1].x;
    verts[i].Verts3D[3].y = -1.0 * verts[i].Verts3D[1].y;
    verts[i].Verts3D[3].z = -1.0 * verts[i].Verts3D[1].z;
    //-----------------------------------------
    for (int j = 0; j < 4; ++j)
    {
        verts[i].Verts3D[j].x += verts[i].Center.x;
        verts[i].Verts3D[j].y += verts[i].Center.y;
        verts[i].Verts3D[j].z += verts[i].Center.z;

        verts[i].Verts3DTruncated[j].x = verts[i].Verts3D[j].x;
        verts[i].Verts3DTruncated[j].y = verts[i].Verts3D[j].y;
        verts[i].Verts3DTruncated[j].z = verts[i].Verts3D[j].z;
    };

    bool gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 0, -1);
    verts[i].ConnectModelSurf[0] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 0, 1);
    verts[i].ConnectModelSurf[1] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 1, -1);
    verts[i].ConnectModelSurf[2] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 1, 1);
    verts[i].ConnectModelSurf[3] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 2, -1);
    verts[i].ConnectModelSurf[4] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 2, 1);
    verts[i].ConnectModelSurf[5] = gh;
}; // FracturesCrossedVertical
template __global__ void
cuDFNsys::FracturesCrossedVertical<double>(cuDFNsys::Fracture<double> *verts,
                                           unsigned long seed, int count,
                                           double model_L);
template __global__ void
cuDFNsys::FracturesCrossedVertical<float>(cuDFNsys::Fracture<float> *verts,
                                          unsigned long seed, int count,
                                          float model_L);

// ====================================================
// NAME:        FracturesBeta50Beta60
// DESCRIPTION: two inclined fractures, with two beta values
// AUTHOR:      Tingchang YIN
// DATE:        20/04/2022
// ====================================================
template <typename T>
__global__ void cuDFNsys::FracturesBeta50Beta60(cuDFNsys::Fracture<T> *verts,
                                                unsigned long seed, int count,
                                                T model_L)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;

    if (i >= 2)
        return;

    curandState state;

    curand_init(seed, i, 0, &state);

    T R_ = model_L * 3.0f;

    verts[i].Radius = R_;
    //printf("%f\n", verts[i].Radius);

    verts[i].NumVertsTruncated = 4;
    //printf("conductivity_powerlaw_exponent: %f, kappa: %f\n", conductivity_powerlaw_exponent, kappa);

    verts[i].Conductivity = (T)(pow(1.0e-3, 3.0) / 12.0);

    if (i == 0)
        verts[i].NormalVec = cuDFNsys::MakeVector3((T)1, (T)0, (T)0);
    else if (i == 1)
        verts[i].NormalVec = cuDFNsys::MakeVector3((T)-1, (T)0, (T)0);

    T R_xy = sqrt(verts[i].NormalVec.x * verts[i].NormalVec.x +
                  verts[i].NormalVec.y * verts[i].NormalVec.y);
    if (i == 0)
        verts[i].NormalVec.z = R_xy / tan(50.0 / 180.0 * M_PI);
    else if (i == 1)
        verts[i].NormalVec.z = R_xy / tan(1.0 / 3.0 * M_PI);
    T norm_f = sqrt(verts[i].NormalVec.x * verts[i].NormalVec.x +
                    verts[i].NormalVec.y * verts[i].NormalVec.y +
                    verts[i].NormalVec.z * verts[i].NormalVec.z);
    verts[i].NormalVec.x /= norm_f;
    verts[i].NormalVec.y /= norm_f;
    verts[i].NormalVec.z /= norm_f;

    verts[i].Center = cuDFNsys::MakeVector3((T)0, (T)0, (T)0);

    //----------
    T *normal_fff = &verts[i].NormalVec.x;
    T *verts_3D_ptr = &verts[i].Verts3D[0].x;
    for (int j = 0; j < 3; ++j)
    {
        if (abs(normal_fff[j]) > 1e-3)
        {
            verts_3D_ptr[(j + 1) % 3] = cuDFNsys::RandomUniform(
                (T)-1.0, (T)1.0, (T)curand_uniform(&state));
            verts_3D_ptr[(j + 2) % 3] = cuDFNsys::RandomUniform(
                (T)-1.0, (T)1.0, (T)curand_uniform(&state));
            verts_3D_ptr[j] =
                -1.0 *
                (verts_3D_ptr[(j + 1) % 3] * normal_fff[(j + 1) % 3] +
                 verts_3D_ptr[(j + 2) % 3] * normal_fff[(j + 2) % 3]) /
                normal_fff[j];
            break;
        }
    }

    T norm_vert1 = sqrt(verts[i].Verts3D[0].x * verts[i].Verts3D[0].x +
                        verts[i].Verts3D[0].y * verts[i].Verts3D[0].y +
                        verts[i].Verts3D[0].z * verts[i].Verts3D[0].z);
    norm_vert1 = R_ / norm_vert1;
    verts[i].Verts3D[0].x *= norm_vert1;
    verts[i].Verts3D[0].y *= norm_vert1;
    verts[i].Verts3D[0].z *= norm_vert1;
    verts[i].Verts3D[2].x = -1.0 * verts[i].Verts3D[0].x;
    verts[i].Verts3D[2].y = -1.0 * verts[i].Verts3D[0].y;
    verts[i].Verts3D[2].z = -1.0 * verts[i].Verts3D[0].z;

    verts[i].Verts3D[1] = cuDFNsys::CrossProductVector3<T>(verts[i].NormalVec,
                                                           verts[i].Verts3D[0]);
    norm_vert1 = sqrt(verts[i].Verts3D[1].x * verts[i].Verts3D[1].x +
                      verts[i].Verts3D[1].y * verts[i].Verts3D[1].y +
                      verts[i].Verts3D[1].z * verts[i].Verts3D[1].z);
    norm_vert1 = R_ / norm_vert1;
    verts[i].Verts3D[1].x *= norm_vert1;
    verts[i].Verts3D[1].y *= norm_vert1;
    verts[i].Verts3D[1].z *= norm_vert1;
    verts[i].Verts3D[3].x = -1.0 * verts[i].Verts3D[1].x;
    verts[i].Verts3D[3].y = -1.0 * verts[i].Verts3D[1].y;
    verts[i].Verts3D[3].z = -1.0 * verts[i].Verts3D[1].z;
    //-----------------------------------------
    for (int j = 0; j < 4; ++j)
    {
        verts[i].Verts3D[j].x += verts[i].Center.x;
        verts[i].Verts3D[j].y += verts[i].Center.y;
        verts[i].Verts3D[j].z += verts[i].Center.z;

        verts[i].Verts3DTruncated[j].x = verts[i].Verts3D[j].x;
        verts[i].Verts3DTruncated[j].y = verts[i].Verts3D[j].y;
        verts[i].Verts3DTruncated[j].z = verts[i].Verts3D[j].z;
    };

    bool gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 0, -1);
    verts[i].ConnectModelSurf[0] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 0, 1);
    verts[i].ConnectModelSurf[1] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 1, -1);
    verts[i].ConnectModelSurf[2] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 1, 1);
    verts[i].ConnectModelSurf[3] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 2, -1);
    verts[i].ConnectModelSurf[4] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 2, 1);
    verts[i].ConnectModelSurf[5] = gh;
}; // FracturesBeta50Beta60
template __global__ void
cuDFNsys::FracturesBeta50Beta60<double>(cuDFNsys::Fracture<double> *verts,
                                        unsigned long seed, int count,
                                        double model_L);
template __global__ void
cuDFNsys::FracturesBeta50Beta60<float>(cuDFNsys::Fracture<float> *verts,
                                       unsigned long seed, int count,
                                       float model_L);

// ====================================================
// NAME:        FracturesIncomplete
// DESCRIPTION: two incomplete fractures
// AUTHOR:      Tingchang YIN
// DATE:        20/04/2022
// ====================================================
template <typename T>
__global__ void cuDFNsys::FracturesIncomplete(cuDFNsys::Fracture<T> *verts,
                                              unsigned long seed, int count,
                                              T model_L)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;

    if (i >= 2)
        return;

    curandState state;

    curand_init(seed, i, 0, &state);

    T R_ = (i == 0 ? 0.55 : 0.60) * model_L;

    verts[i].Radius = R_;
    //printf("%f\n", verts[i].Radius);

    verts[i].NumVertsTruncated = 4;
    //printf("conductivity_powerlaw_exponent: %f, kappa: %f\n", conductivity_powerlaw_exponent, kappa);

    verts[i].Conductivity = (T)(pow(1.0e-3, 3.0) / 12.0);

    if (i == 0)
        verts[i].NormalVec = cuDFNsys::MakeVector3((T)1, (T)0, (T)0);
    else if (i == 1)
        verts[i].NormalVec = cuDFNsys::MakeVector3((T)0, (T)1, (T)0);

    verts[i].Center = cuDFNsys::MakeVector3((T)0, (T)0, (T)0);

    //----------
    if (i == 0)
    {
        verts[i].Verts3D[0].y = 0.5 * model_L;
        verts[i].Verts3D[0].x = 0;
        verts[i].Verts3D[0].z = 0;
    }
    else if (i == 1)
    {
        verts[i].Verts3D[0].x = 0.5 * model_L;
        verts[i].Verts3D[0].y = 0;
        verts[i].Verts3D[0].z = 0;
    }

    T norm_vert1 = sqrt(verts[i].Verts3D[0].x * verts[i].Verts3D[0].x +
                        verts[i].Verts3D[0].y * verts[i].Verts3D[0].y +
                        verts[i].Verts3D[0].z * verts[i].Verts3D[0].z);
    norm_vert1 = R_ / norm_vert1;
    verts[i].Verts3D[0].x *= norm_vert1;
    verts[i].Verts3D[0].y *= norm_vert1;
    verts[i].Verts3D[0].z *= norm_vert1;
    verts[i].Verts3D[2].x = -1.0 * verts[i].Verts3D[0].x;
    verts[i].Verts3D[2].y = -1.0 * verts[i].Verts3D[0].y;
    verts[i].Verts3D[2].z = -1.0 * verts[i].Verts3D[0].z;

    verts[i].Verts3D[1] = cuDFNsys::CrossProductVector3<T>(verts[i].NormalVec,
                                                           verts[i].Verts3D[0]);
    norm_vert1 = sqrt(verts[i].Verts3D[1].x * verts[i].Verts3D[1].x +
                      verts[i].Verts3D[1].y * verts[i].Verts3D[1].y +
                      verts[i].Verts3D[1].z * verts[i].Verts3D[1].z);
    norm_vert1 = R_ / norm_vert1;
    verts[i].Verts3D[1].x *= norm_vert1;
    verts[i].Verts3D[1].y *= norm_vert1;
    verts[i].Verts3D[1].z *= norm_vert1;
    verts[i].Verts3D[3].x = -1.0 * verts[i].Verts3D[1].x;
    verts[i].Verts3D[3].y = -1.0 * verts[i].Verts3D[1].y;
    verts[i].Verts3D[3].z = -1.0 * verts[i].Verts3D[1].z;
    //-----------------------------------------
    for (int j = 0; j < 4; ++j)
    {
        verts[i].Verts3D[j].x += verts[i].Center.x;
        verts[i].Verts3D[j].y += verts[i].Center.y;
        verts[i].Verts3D[j].z += verts[i].Center.z;

        verts[i].Verts3DTruncated[j].x = verts[i].Verts3D[j].x;
        verts[i].Verts3DTruncated[j].y = verts[i].Verts3D[j].y;
        verts[i].Verts3DTruncated[j].z = verts[i].Verts3D[j].z;
    };

    bool gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 0, -1);
    verts[i].ConnectModelSurf[0] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 0, 1);
    verts[i].ConnectModelSurf[1] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 1, -1);
    verts[i].ConnectModelSurf[2] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 1, 1);
    verts[i].ConnectModelSurf[3] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 2, -1);
    verts[i].ConnectModelSurf[4] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 2, 1);
    verts[i].ConnectModelSurf[5] = gh;
}; // FracturesIncomplete
template __global__ void
cuDFNsys::FracturesIncomplete<double>(cuDFNsys::Fracture<double> *verts,
                                      unsigned long seed, int count,
                                      double model_L);
template __global__ void
cuDFNsys::FracturesIncomplete<float>(cuDFNsys::Fracture<float> *verts,
                                     unsigned long seed, int count,
                                     float model_L);

// ====================================================
// NAME:        Fractures2DLike
// DESCRIPTION: 2D-like fractures
// AUTHOR:      Tingchang YIN
// DATE:        20/04/2022
// ====================================================
template <typename T>
__global__ void cuDFNsys::Fractures2DLike(cuDFNsys::Fracture<T> *verts,
                                          unsigned long seed, int count,
                                          T model_L, T alpha, T minR, T maxR)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;

    if (i > count - 1)
        return;

    curandState state;

    curand_init(seed, i, 0, &state);

    T R_ = 0;

    R_ = cuDFNsys::RandomPowerlaw((T)minR, (T)maxR, (T)alpha,
                                  (T)curand_uniform(&state));

    verts[i].Radius = R_;
    //printf("%f\n", verts[i].Radius);

    verts[i].NumVertsTruncated = 4;
    //printf("conductivity_powerlaw_exponent: %f, kappa: %f\n", conductivity_powerlaw_exponent, kappa);

    verts[i].Conductivity = (T)(pow(1.0e-3, 3.0) / 12.0);

    verts[i].Center.x = cuDFNsys::RandomUniform(
        (T)(-model_L * 0.5), (T)(model_L * 0.5), (T)(curand_uniform(&state)));
    verts[i].Center.y = 0;
    verts[i].Center.z = cuDFNsys::RandomUniform(
        (T)(-model_L * 0.5), (T)(model_L * 0.5), (T)(curand_uniform(&state)));

    verts[i].NormalVec = cuDFNsys::MakeVector3(
        (T)cuDFNsys::RandomUniform((T)-1.0, (T)1.0, (T)curand_uniform(&state)),
        (T)0,
        (T)cuDFNsys::RandomUniform((T)-1.0, (T)1.0, (T)curand_uniform(&state)));
    T norm_f = sqrt(verts[i].NormalVec.x * verts[i].NormalVec.x +
                    verts[i].NormalVec.y * verts[i].NormalVec.y +
                    verts[i].NormalVec.z * verts[i].NormalVec.z);
    verts[i].NormalVec.x /= norm_f;
    verts[i].NormalVec.y /= norm_f;
    verts[i].NormalVec.z /= norm_f;

    //--------------
    cuDFNsys::Vector3<T> rotate_axis = cuDFNsys::MakeVector3(
        -verts[i].NormalVec.y, verts[i].NormalVec.x, (T)0);
    T rotate_axis_norm =
        sqrt(rotate_axis.x * rotate_axis.x + rotate_axis.y * rotate_axis.y +
             rotate_axis.z * rotate_axis.z);
    rotate_axis.x /= rotate_axis_norm;
    rotate_axis.y /= rotate_axis_norm;
    rotate_axis.z /= rotate_axis_norm;

    cuDFNsys::Quaternion<T> qua;
    qua = qua.DescribeRotation(verts[i].NormalVec, 45.0 / 180.0 * M_PI);
    verts[i].Verts3D[0] = qua.Rotate(rotate_axis);

    //--------------
    T norm_vert1 = sqrt(verts[i].Verts3D[0].x * verts[i].Verts3D[0].x +
                        verts[i].Verts3D[0].y * verts[i].Verts3D[0].y +
                        verts[i].Verts3D[0].z * verts[i].Verts3D[0].z);
    norm_vert1 = R_ / norm_vert1;
    verts[i].Verts3D[0].x *= norm_vert1;
    verts[i].Verts3D[0].y *= norm_vert1;
    verts[i].Verts3D[0].z *= norm_vert1;
    verts[i].Verts3D[2].x = -1.0 * verts[i].Verts3D[0].x;
    verts[i].Verts3D[2].y = -1.0 * verts[i].Verts3D[0].y;
    verts[i].Verts3D[2].z = -1.0 * verts[i].Verts3D[0].z;

    verts[i].Verts3D[1] = cuDFNsys::CrossProductVector3<T>(verts[i].NormalVec,
                                                           verts[i].Verts3D[0]);
    norm_vert1 = sqrt(verts[i].Verts3D[1].x * verts[i].Verts3D[1].x +
                      verts[i].Verts3D[1].y * verts[i].Verts3D[1].y +
                      verts[i].Verts3D[1].z * verts[i].Verts3D[1].z);
    norm_vert1 = R_ / norm_vert1;
    verts[i].Verts3D[1].x *= norm_vert1;
    verts[i].Verts3D[1].y *= norm_vert1;
    verts[i].Verts3D[1].z *= norm_vert1;
    verts[i].Verts3D[3].x = -1.0 * verts[i].Verts3D[1].x;
    verts[i].Verts3D[3].y = -1.0 * verts[i].Verts3D[1].y;
    verts[i].Verts3D[3].z = -1.0 * verts[i].Verts3D[1].z;
    //-----------------------------------------
    for (int j = 0; j < 4; ++j)
    {
        verts[i].Verts3D[j].x += verts[i].Center.x;
        verts[i].Verts3D[j].y += verts[i].Center.y;
        verts[i].Verts3D[j].z += verts[i].Center.z;

        verts[i].Verts3D[j].y =
            (verts[i].Verts3D[j].y / abs(verts[i].Verts3D[j].y)) * 0.5 *
            model_L;
        verts[i].Verts3DTruncated[j].x = verts[i].Verts3D[j].x;
        verts[i].Verts3DTruncated[j].y = verts[i].Verts3D[j].y;
        verts[i].Verts3DTruncated[j].z = verts[i].Verts3D[j].z;
    };

    bool gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 0, -1);
    verts[i].ConnectModelSurf[0] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 0, 1);
    verts[i].ConnectModelSurf[1] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 1, -1);
    verts[i].ConnectModelSurf[2] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 1, 1);
    verts[i].ConnectModelSurf[3] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 2, -1);
    verts[i].ConnectModelSurf[4] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 2, 1);
    verts[i].ConnectModelSurf[5] = gh;
}; // Fractures2DLike
template __global__ void
cuDFNsys::Fractures2DLike<double>(cuDFNsys::Fracture<double> *verts,
                                  unsigned long seed, int count, double model_L,
                                  double alpha, double minR, double maxR);
template __global__ void
cuDFNsys::Fractures2DLike<float>(cuDFNsys::Fracture<float> *verts,
                                 unsigned long seed, int count, float model_L,
                                 float alpha, float minR, float maxR);

// ====================================================
// NAME:        FracturesFour
// DESCRIPTION: Four fractures to verify particle tracking
// AUTHOR:      Tingchang YIN
// DATE:        01/09/2022
// ====================================================
template <typename T>
__global__ void cuDFNsys::FracturesFour(cuDFNsys::Fracture<T> *verts,
                                        unsigned long seed, int count,
                                        T model_L)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;

    if (i >= 4)
        return;

    curandState state;

    curand_init(seed, i, 0, &state);

    T R_ = sqrt(model_L * model_L / 2.0f);

    verts[i].Radius = R_;
    //printf("%f\n", verts[i].Radius);

    verts[i].NumVertsTruncated = 4;
    //printf("conductivity_powerlaw_exponent: %f, kappa: %f\n", conductivity_powerlaw_exponent, kappa);

    verts[i].Conductivity = (T)(pow(1.0e-3, 3.0) / 12.0);

    if (i == 0)
        verts[i].Center = cuDFNsys::MakeVector3((T)0., (T)0., (T)0.); //
    else if (i == 1)
        verts[i].Center =
            cuDFNsys::MakeVector3((T)0., (T)0., (T)0.25 * model_L);
    else if (i == 2)
        verts[i].Center = cuDFNsys::MakeVector3((T)0.25 * model_L, (T)0.,
                                                (T)(-0.25) * model_L);
    else if (i == 3)
        verts[i].Center = cuDFNsys::MakeVector3((T)(-0.25) * model_L, (T)0.,
                                                (T)(-0.25) * model_L);

    if (i == 0)
        verts[i].NormalVec = cuDFNsys::MakeVector3((T)0., (T)0., (T)1.); //
    else if (i == 1)
        verts[i].NormalVec = cuDFNsys::MakeVector3(-(T)1., (T)0., (T)0.);
    else if (i == 2)
        verts[i].NormalVec = cuDFNsys::MakeVector3(-(T)1., (T)0., (T)0.);
    else if (i == 3)
        verts[i].NormalVec = cuDFNsys::MakeVector3(-(T)1., (T)0., (T)0.);

    if (i == 0)
        verts[i].Verts3D[0] =
            cuDFNsys::MakeVector3((T)0.5 * model_L, (T)0.5 * model_L, (T)0.); //
    else
        verts[i].Verts3D[0] =
            cuDFNsys::MakeVector3((T)0., (T)0.5 * model_L, (T)(0.5) * model_L);

    T norm_vert1 = sqrt(verts[i].Verts3D[0].x * verts[i].Verts3D[0].x +
                        verts[i].Verts3D[0].y * verts[i].Verts3D[0].y +
                        verts[i].Verts3D[0].z * verts[i].Verts3D[0].z);
    norm_vert1 = R_ / norm_vert1;
    verts[i].Verts3D[0].x *= norm_vert1;
    verts[i].Verts3D[0].y *= norm_vert1;
    verts[i].Verts3D[0].z *= norm_vert1;
    verts[i].Verts3D[2].x = -1.0 * verts[i].Verts3D[0].x;
    verts[i].Verts3D[2].y = -1.0 * verts[i].Verts3D[0].y;
    verts[i].Verts3D[2].z = -1.0 * verts[i].Verts3D[0].z;

    verts[i].Verts3D[1] = cuDFNsys::CrossProductVector3<T>(verts[i].NormalVec,
                                                           verts[i].Verts3D[0]);
    norm_vert1 = sqrt(verts[i].Verts3D[1].x * verts[i].Verts3D[1].x +
                      verts[i].Verts3D[1].y * verts[i].Verts3D[1].y +
                      verts[i].Verts3D[1].z * verts[i].Verts3D[1].z);
    norm_vert1 = R_ / norm_vert1;
    verts[i].Verts3D[1].x *= norm_vert1;
    verts[i].Verts3D[1].y *= norm_vert1;
    verts[i].Verts3D[1].z *= norm_vert1;
    verts[i].Verts3D[3].x = -1.0 * verts[i].Verts3D[1].x;
    verts[i].Verts3D[3].y = -1.0 * verts[i].Verts3D[1].y;
    verts[i].Verts3D[3].z = -1.0 * verts[i].Verts3D[1].z;

    //-----------------------------------------
    for (int j = 0; j < 4; ++j)
    {
        verts[i].Verts3D[j].x += verts[i].Center.x;
        verts[i].Verts3D[j].y += verts[i].Center.y;
        verts[i].Verts3D[j].z += verts[i].Center.z;

        verts[i].Verts3DTruncated[j].x = verts[i].Verts3D[j].x;
        verts[i].Verts3DTruncated[j].y = verts[i].Verts3D[j].y;
        verts[i].Verts3DTruncated[j].z = verts[i].Verts3D[j].z;
    };

    bool gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 0, -1);
    verts[i].ConnectModelSurf[0] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 0, 1);
    verts[i].ConnectModelSurf[1] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 1, -1);
    verts[i].ConnectModelSurf[2] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 1, 1);
    verts[i].ConnectModelSurf[3] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 2, -1);
    verts[i].ConnectModelSurf[4] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 2, 1);
    verts[i].ConnectModelSurf[5] = gh;
};
template __global__ void
cuDFNsys::FracturesFour<double>(cuDFNsys::Fracture<double> *verts,
                                unsigned long seed, int count, double model_L);
template __global__ void
cuDFNsys::FracturesFour<float>(cuDFNsys::Fracture<float> *verts,
                               unsigned long seed, int count, float model_L);

// ====================================================
// NAME:        FracturesChangeDomainSize
// DESCRIPTION: Change Domain size
// AUTHOR:      Tingchang YIN
// DATE:        13/09/2022
// ====================================================
template <typename T>
__global__ void
cuDFNsys::FracturesChangeDomainSize(cuDFNsys::Fracture<T> *verts, int count,
                                    T model_L)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;

    if (i > count - 1)
        return;

    verts[i].NumVertsTruncated = 4;
    for (int j = 0; j < 4; ++j)
        verts[i].Verts3DTruncated[j].x = verts[i].Verts3D[j].x,
        verts[i].Verts3DTruncated[j].y = verts[i].Verts3D[j].y,
        verts[i].Verts3DTruncated[j].z = verts[i].Verts3D[j].z;

    bool gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 0, -1);
    verts[i].ConnectModelSurf[0] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 0, 1);
    verts[i].ConnectModelSurf[1] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 1, -1);
    verts[i].ConnectModelSurf[2] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 1, 1);
    verts[i].ConnectModelSurf[3] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 2, -1);
    verts[i].ConnectModelSurf[4] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 2, 1);
    verts[i].ConnectModelSurf[5] = gh;
}; // FracturesChangeDomainSize
template __global__ void
cuDFNsys::FracturesChangeDomainSize<double>(cuDFNsys::Fracture<double> *verts,
                                            int count, double model_L);
template __global__ void
cuDFNsys::FracturesChangeDomainSize<float>(cuDFNsys::Fracture<float> *verts,
                                           int count, float model_L);

// ====================================================
// NAME:        FractureTwoIntersectOrNot
// DESCRIPTION: two fractures intersect or not
// AUTHOR:      Tingchang YIN
// DATE:        13/06/2023
// ====================================================
template <typename T>
__global__ void
cuDFNsys::FractureTwoIntersectOrNot(cuDFNsys::Fracture<T> *verts,
                                    unsigned long seed, int count, T model_L,
                                    bool IfIntersect)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;

    if (i > 1)
        return;

    curandState state;
    // printf("i = %d\n", i);
    curand_init(seed, i, 0, &state);

    T R_ = model_L * 2;

    verts[i].Radius = R_;
    //printf("%f\n", verts[i].Radius);

    verts[i].NumVertsTruncated = 4;
    //printf("conductivity_powerlaw_exponent: %f, kappa: %f\n", conductivity_powerlaw_exponent, kappa);

    if (i == 0)
        verts[i].Conductivity = (T)(pow(1.0e-3, 3.0) / 12.0);
    else
        verts[i].Conductivity = (T)(pow(2.0e-3, 3.0) / 12.0);

    if (i == 0)
        verts[i].NormalVec = cuDFNsys::MakeVector3(
            (T)0., -(T)(1.0 / tan(20.f / 180.f * M_PI)), (T)1.);
    else
        verts[i].NormalVec = cuDFNsys::MakeVector3(
            (T)0., (T)(1.0 / tan(20.f / 180.f * M_PI)), (T)1.);

    T norm_NormalVec = sqrt(verts[i].NormalVec.x * verts[i].NormalVec.x +
                            verts[i].NormalVec.y * verts[i].NormalVec.y +
                            verts[i].NormalVec.z * verts[i].NormalVec.z);
    verts[i].NormalVec.x /= norm_NormalVec;
    verts[i].NormalVec.y /= norm_NormalVec;
    verts[i].NormalVec.z /= norm_NormalVec;

    if (IfIntersect)
        verts[i].Center = cuDFNsys::MakeVector3((T)0., (T)0., (T)0.);
    else
    {
        if (i == 0)
            verts[i].Center = cuDFNsys::MakeVector3(
                (T)0., -(T)(0.55f * model_L * tan(20.f / 180.f * M_PI)), (T)0.);
        else
            verts[i].Center = cuDFNsys::MakeVector3(
                (T)0., (T)(0.55f * model_L * tan(20.f / 180.f * M_PI)), (T)0.);
    }

    cuDFNsys::Vector1<T> *normal_fff = &verts[i].NormalVec.x;
    cuDFNsys::Vector1<T> *verts_3D_ptr = &verts[i].Verts3D[0].x;
    for (int j = 0; j < 3; ++j)
    {
        if (abs(normal_fff[j]) > 1e-3)
        {
            verts_3D_ptr[(j + 1) % 3] = cuDFNsys::RandomUniform(
                (cuDFNsys::Vector1<T>)-1.0, (cuDFNsys::Vector1<T>)1.0,
                (cuDFNsys::Vector1<T>)curand_uniform(&state));
            verts_3D_ptr[(j + 2) % 3] = cuDFNsys::RandomUniform(
                (cuDFNsys::Vector1<T>)-1.0, (cuDFNsys::Vector1<T>)1.0,
                (cuDFNsys::Vector1<T>)curand_uniform(&state));
            verts_3D_ptr[j] =
                -1.0 *
                (verts_3D_ptr[(j + 1) % 3] * normal_fff[(j + 1) % 3] +
                 verts_3D_ptr[(j + 2) % 3] * normal_fff[(j + 2) % 3]) /
                normal_fff[j];
            break;
        }
    }

    cuDFNsys::Vector1<T> norm_vert1 =
        sqrt(verts[i].Verts3D[0].x * verts[i].Verts3D[0].x +
             verts[i].Verts3D[0].y * verts[i].Verts3D[0].y +
             verts[i].Verts3D[0].z * verts[i].Verts3D[0].z);
    norm_vert1 = R_ / norm_vert1;
    verts[i].Verts3D[0].x *= norm_vert1;
    verts[i].Verts3D[0].y *= norm_vert1;
    verts[i].Verts3D[0].z *= norm_vert1;
    verts[i].Verts3D[2].x = -1.0 * verts[i].Verts3D[0].x;
    verts[i].Verts3D[2].y = -1.0 * verts[i].Verts3D[0].y;
    verts[i].Verts3D[2].z = -1.0 * verts[i].Verts3D[0].z;

    verts[i].Verts3D[1] = cuDFNsys::CrossProductVector3<T>(verts[i].NormalVec,
                                                           verts[i].Verts3D[0]);
    norm_vert1 = sqrt(verts[i].Verts3D[1].x * verts[i].Verts3D[1].x +
                      verts[i].Verts3D[1].y * verts[i].Verts3D[1].y +
                      verts[i].Verts3D[1].z * verts[i].Verts3D[1].z);
    norm_vert1 = R_ / norm_vert1;
    verts[i].Verts3D[1].x *= norm_vert1;
    verts[i].Verts3D[1].y *= norm_vert1;
    verts[i].Verts3D[1].z *= norm_vert1;
    verts[i].Verts3D[3].x = -1.0 * verts[i].Verts3D[1].x;
    verts[i].Verts3D[3].y = -1.0 * verts[i].Verts3D[1].y;
    verts[i].Verts3D[3].z = -1.0 * verts[i].Verts3D[1].z;
    //-----------------------------------------
    for (int j = 0; j < 4; ++j)
    {
        verts[i].Verts3D[j].x += verts[i].Center.x;
        verts[i].Verts3D[j].y += verts[i].Center.y;
        verts[i].Verts3D[j].z += verts[i].Center.z;

        verts[i].Verts3DTruncated[j].x = verts[i].Verts3D[j].x;
        verts[i].Verts3DTruncated[j].y = verts[i].Verts3D[j].y;
        verts[i].Verts3DTruncated[j].z = verts[i].Verts3D[j].z;
    };

    bool gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 0, -1);
    verts[i].ConnectModelSurf[0] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 0, 1);
    verts[i].ConnectModelSurf[1] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 1, -1);
    verts[i].ConnectModelSurf[2] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 1, 1);
    verts[i].ConnectModelSurf[3] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 2, -1);
    verts[i].ConnectModelSurf[4] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 2, 1);
    verts[i].ConnectModelSurf[5] = gh;
};
template __global__ void
cuDFNsys::FractureTwoIntersectOrNot<double>(cuDFNsys::Fracture<double> *verts,
                                            unsigned long seed, int count,
                                            double model_L, bool IfIntersect);
template __global__ void
cuDFNsys::FractureTwoIntersectOrNot<float>(cuDFNsys::Fracture<float> *verts,
                                           unsigned long seed, int count,
                                           float model_L, bool IfIntersect);

template <typename T>
__global__ void cuDFNsys::FracturesParallel(cuDFNsys::Fracture<T> *verts,
                                            int count, unsigned long seed,
                                            T model_L)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;

    if (i >= count)
        return;

    curandState state;
    // printf("i = %d\n", i);
    curand_init(seed, i, 0, &state);

    T R_ = model_L * 2;

    verts[i].Radius = R_;
    //printf("%f\n", verts[i].Radius);

    verts[i].NumVertsTruncated = 4;
    //printf("conductivity_powerlaw_exponent: %f, kappa: %f\n", conductivity_powerlaw_exponent, kappa);

    verts[i].Conductivity = cuDFNsys::RandomUniform(
        (T)1e-10, (T)1e-6, (cuDFNsys::Vector1<T>)curand_uniform(&state));

    // verts[i].NormalVec = cuDFNsys::MakeVector3((T)0., (T)1, (T)0.);
    verts[i].NormalVec = cuDFNsys::MakeVector3(
        (T)0.,
        (T)1, //cuDFNsys::RandomUniform((T)-1, (T)1, (cuDFNsys::Vector1<T>)curand_uniform(&state)),
        cuDFNsys::RandomUniform((T)-1, (T)1,
                                (cuDFNsys::Vector1<T>)curand_uniform(&state)));

    if (verts[i].NormalVec.z < 0)
        verts[i].NormalVec = cuDFNsys::MakeVector3((T)0, -verts[i].NormalVec.y,
                                                   -verts[i].NormalVec.z);

    T norm_NormalVec = sqrt(verts[i].NormalVec.x * verts[i].NormalVec.x +
                            verts[i].NormalVec.y * verts[i].NormalVec.y +
                            verts[i].NormalVec.z * verts[i].NormalVec.z);
    verts[i].NormalVec.x /= norm_NormalVec;
    verts[i].NormalVec.y /= norm_NormalVec;
    verts[i].NormalVec.z /= norm_NormalVec;

    // rintf("%f, %f, %f \n", verts[i].NormalVec.x, verts[i].NormalVec.y, verts[i].NormalVec.z);

    T spacing = model_L / (count + 1);

    verts[i].Center = cuDFNsys::MakeVector3(
        (T)0., (T)(-0.5 * model_L + (i + 1) * spacing), (T)0.);

    cuDFNsys::Vector1<T> *normal_fff = &verts[i].NormalVec.x;
    cuDFNsys::Vector1<T> *verts_3D_ptr = &verts[i].Verts3D[0].x;
    for (int j = 0; j < 3; ++j)
    {
        if (abs(normal_fff[j]) > 1e-3)
        {
            verts_3D_ptr[(j + 1) % 3] = cuDFNsys::RandomUniform(
                (cuDFNsys::Vector1<T>)-1.0, (cuDFNsys::Vector1<T>)1.0,
                (cuDFNsys::Vector1<T>)curand_uniform(&state));
            verts_3D_ptr[(j + 2) % 3] = cuDFNsys::RandomUniform(
                (cuDFNsys::Vector1<T>)-1.0, (cuDFNsys::Vector1<T>)1.0,
                (cuDFNsys::Vector1<T>)curand_uniform(&state));
            verts_3D_ptr[j] =
                -1.0 *
                (verts_3D_ptr[(j + 1) % 3] * normal_fff[(j + 1) % 3] +
                 verts_3D_ptr[(j + 2) % 3] * normal_fff[(j + 2) % 3]) /
                normal_fff[j];
            break;
        }
    }

    cuDFNsys::Vector1<T> norm_vert1 =
        sqrt(verts[i].Verts3D[0].x * verts[i].Verts3D[0].x +
             verts[i].Verts3D[0].y * verts[i].Verts3D[0].y +
             verts[i].Verts3D[0].z * verts[i].Verts3D[0].z);
    norm_vert1 = R_ / norm_vert1;
    verts[i].Verts3D[0].x *= norm_vert1;
    verts[i].Verts3D[0].y *= norm_vert1;
    verts[i].Verts3D[0].z *= norm_vert1;
    verts[i].Verts3D[2].x = -1.0 * verts[i].Verts3D[0].x;
    verts[i].Verts3D[2].y = -1.0 * verts[i].Verts3D[0].y;
    verts[i].Verts3D[2].z = -1.0 * verts[i].Verts3D[0].z;

    verts[i].Verts3D[1] = cuDFNsys::CrossProductVector3<T>(verts[i].NormalVec,
                                                           verts[i].Verts3D[0]);
    norm_vert1 = sqrt(verts[i].Verts3D[1].x * verts[i].Verts3D[1].x +
                      verts[i].Verts3D[1].y * verts[i].Verts3D[1].y +
                      verts[i].Verts3D[1].z * verts[i].Verts3D[1].z);
    norm_vert1 = R_ / norm_vert1;
    verts[i].Verts3D[1].x *= norm_vert1;
    verts[i].Verts3D[1].y *= norm_vert1;
    verts[i].Verts3D[1].z *= norm_vert1;
    verts[i].Verts3D[3].x = -1.0 * verts[i].Verts3D[1].x;
    verts[i].Verts3D[3].y = -1.0 * verts[i].Verts3D[1].y;
    verts[i].Verts3D[3].z = -1.0 * verts[i].Verts3D[1].z;
    //-----------------------------------------
    for (int j = 0; j < 4; ++j)
    {
        verts[i].Verts3D[j].x += verts[i].Center.x;
        verts[i].Verts3D[j].y += verts[i].Center.y;
        verts[i].Verts3D[j].z += verts[i].Center.z;

        verts[i].Verts3DTruncated[j].x = verts[i].Verts3D[j].x;
        verts[i].Verts3DTruncated[j].y = verts[i].Verts3D[j].y;
        verts[i].Verts3DTruncated[j].z = verts[i].Verts3D[j].z;
    };

    bool gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 0, -1);
    verts[i].ConnectModelSurf[0] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 0, 1);
    verts[i].ConnectModelSurf[1] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 1, -1);
    verts[i].ConnectModelSurf[2] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 1, 1);
    verts[i].ConnectModelSurf[3] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 2, -1);
    verts[i].ConnectModelSurf[4] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 2, 1);
    verts[i].ConnectModelSurf[5] = gh;
}; // FracturesParallel
template __global__ void
cuDFNsys::FracturesParallel<double>(cuDFNsys::Fracture<double> *verts,
                                    int count, unsigned long seed,
                                    double model_L);
template __global__ void
cuDFNsys::FracturesParallel<float>(cuDFNsys::Fracture<float> *verts, int count,
                                   unsigned long seed, float model_L);