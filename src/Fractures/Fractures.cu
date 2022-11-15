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
__global__ void cuDFNsys::Fractures(cuDFNsys::Fracture<T> *verts,
                                    unsigned long seed,
                                    int count,
                                    cuDFNsys::Vector1<T> model_L,
                                    uint ModeSizeDistri,                 // 1 = power law; 2 = lognormal; 3 = uniform; 4 = monosize
                                    cuDFNsys::Vector4<T> ParaSizeDistri, // when mode = 1, ;
                                    cuDFNsys::Vector1<T> kappa,
                                    cuDFNsys::Vector1<T> conductivity_powerlaw_exponent)
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
                                      (T)ParaSizeDistri.x, (T)curand_uniform(&state));
    else if (ModeSizeDistri == 1)
        R_ = cuDFNsys::RandomLognormal((T)ParaSizeDistri.x,
                                       (T)ParaSizeDistri.y,
                                       (T)ParaSizeDistri.z,
                                       (T)ParaSizeDistri.w, (T)curand_uniform(&state));
    else if (ModeSizeDistri == 2)
        R_ = cuDFNsys::RandomUniform((T)ParaSizeDistri.x,
                                     (T)ParaSizeDistri.y, (T)curand_uniform(&state));
    else if (ModeSizeDistri == 3)
        R_ = ParaSizeDistri.x;

    verts[i].Radius = R_;
    //printf("%f\n", verts[i].Radius);

    verts[i].NumVertsTruncated = 4;
    //printf("conductivity_powerlaw_exponent: %f, kappa: %f\n", conductivity_powerlaw_exponent, kappa);

    if (conductivity_powerlaw_exponent == 0)
        verts[i].Conductivity = (T)(pow(1.0e-3, 3.0) / 12.0);
    else
        verts[i].Conductivity = (1.0e-11) * pow(R_, 3.0 * conductivity_powerlaw_exponent);

    verts[i].Center.x = cuDFNsys::RandomUniform((T)(-model_L * 0.5),
                                                (T)(model_L * 0.5),
                                                (T)(curand_uniform(&state)));
    verts[i].Center.y = cuDFNsys::RandomUniform((T)(-model_L * 0.5),
                                                (T)(model_L * 0.5),
                                                (T)(curand_uniform(&state)));
    verts[i].Center.z = cuDFNsys::RandomUniform((T)(-model_L * 0.5),
                                                (T)(model_L * 0.5),
                                                (T)(curand_uniform(&state)));

    verts[i].NormalVec = cuDFNsys::MakeVector3((T)cuDFNsys::RandomUniform((cuDFNsys::Vector1<T>)-1.0, (cuDFNsys::Vector1<T>)1.0, (cuDFNsys::Vector1<T>)curand_uniform(&state)),
                                               (T)cuDFNsys::RandomUniform((cuDFNsys::Vector1<T>)-1.0, (cuDFNsys::Vector1<T>)1.0, (cuDFNsys::Vector1<T>)curand_uniform(&state)),
                                               (T)0);
    cuDFNsys::Vector1<T> R_xy = sqrt(verts[i].NormalVec.x * verts[i].NormalVec.x +
                                     verts[i].NormalVec.y * verts[i].NormalVec.y);
    verts[i].NormalVec.z = R_xy / tan(cuDFNsys::RandomFisher((T)curand_uniform(&state),
                                                             (T)kappa));
    cuDFNsys::Vector1<T> norm_f = sqrt(verts[i].NormalVec.x * verts[i].NormalVec.x +
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
            verts_3D_ptr[(j + 1) % 3] = cuDFNsys::RandomUniform((cuDFNsys::Vector1<T>)-1.0, (cuDFNsys::Vector1<T>)1.0, (cuDFNsys::Vector1<T>)curand_uniform(&state));
            verts_3D_ptr[(j + 2) % 3] = cuDFNsys::RandomUniform((cuDFNsys::Vector1<T>)-1.0, (cuDFNsys::Vector1<T>)1.0, (cuDFNsys::Vector1<T>)curand_uniform(&state));
            verts_3D_ptr[j] = -1.0 * (verts_3D_ptr[(j + 1) % 3] * normal_fff[(j + 1) % 3] + verts_3D_ptr[(j + 2) % 3] * normal_fff[(j + 2) % 3]) / normal_fff[j];
            break;
        }
    }

    cuDFNsys::Vector1<T> norm_vert1 = sqrt(verts[i].Verts3D[0].x * verts[i].Verts3D[0].x +
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

    bool gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 0, 1);
    verts[i].ConnectModelSurf[0] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 0, -1);
    verts[i].ConnectModelSurf[1] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 1, 1);
    verts[i].ConnectModelSurf[2] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 1, -1);
    verts[i].ConnectModelSurf[3] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 2, 1);
    verts[i].ConnectModelSurf[4] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 2, -1);
    verts[i].ConnectModelSurf[5] = gh;
}; // Fractures
template __global__ void cuDFNsys::Fractures<double>(cuDFNsys::Fracture<double> *verts,
                                                     unsigned long seed,
                                                     int count,
                                                     cuDFNsys::Vector1<double> model_L,
                                                     uint ModeSizeDistri,                      // 1 = power law; 2 = lognormal; 3 = uniform; 4 = monosize
                                                     cuDFNsys::Vector4<double> ParaSizeDistri, // when mode = 1, ;
                                                     cuDFNsys::Vector1<double> kappa,
                                                     cuDFNsys::Vector1<double> conductivity_powerlaw_exponent);
template __global__ void cuDFNsys::Fractures<float>(cuDFNsys::Fracture<float> *verts,
                                                    unsigned long seed,
                                                    int count,
                                                    cuDFNsys::Vector1<float> model_L,
                                                    uint ModeSizeDistri,                     // 1 = power law; 2 = lognormal; 3 = uniform; 4 = monosize
                                                    cuDFNsys::Vector4<float> ParaSizeDistri, // when mode = 1, ;
                                                    cuDFNsys::Vector1<float> kappa,
                                                    cuDFNsys::Vector1<float> conductivity_powerlaw_exponent);

// ====================================================
// NAME:        Fractures
// DESCRIPTION: Fractures in a DFN
// AUTHOR:      Tingchang YIN
// DATE:        04/04/2022
// ====================================================
template <typename T>
cuDFNsys::FracturesCPU<T>::FracturesCPU(thrust::host_vector<cuDFNsys::Fracture<T>> &verts,
                                        unsigned long seed,
                                        int count,
                                        T model_L,
                                        uint ModeSizeDistri,                 // 0 = power law; 1 = lognormal; 2 = uniform; 3 = monosize
                                        cuDFNsys::Vector4<T> ParaSizeDistri, // when mode = 1, ;
                                        T kappa,
                                        T conductivity_powerlaw_exponent, uint Nproc)
{
    srand(seed);

    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (int i = 0; i < count; ++i)
    {
        T R_ = 0;

        if (ModeSizeDistri == 0)
            R_ = cuDFNsys::RandomPowerlaw((T)ParaSizeDistri.y, (T)ParaSizeDistri.z,
                                          (T)ParaSizeDistri.x, (T)((T)rand() / (T)RAND_MAX));
        else if (ModeSizeDistri == 1)
            R_ = cuDFNsys::RandomLognormal((T)ParaSizeDistri.x,
                                           (T)ParaSizeDistri.y,
                                           (T)ParaSizeDistri.z,
                                           (T)ParaSizeDistri.w, (T)((T)rand() / (T)RAND_MAX));
        else if (ModeSizeDistri == 2)
            R_ = cuDFNsys::RandomUniform((T)ParaSizeDistri.x,
                                         (T)ParaSizeDistri.y, (T)((T)rand() / (T)RAND_MAX));
        else if (ModeSizeDistri == 3)
            R_ = ParaSizeDistri.x;

        verts[i].Radius = R_;
        //printf("%f\n", verts[i].Radius);

        verts[i].NumVertsTruncated = 4;
        //printf("conductivity_powerlaw_exponent: %f, kappa: %f\n", conductivity_powerlaw_exponent, kappa);

        if (conductivity_powerlaw_exponent == 0)
            verts[i].Conductivity = (T)(pow(1.0e-3, 3.0) / 12.0);
        else
            verts[i].Conductivity = (1.0e-11) * pow(R_, 3.0 * conductivity_powerlaw_exponent);

        verts[i].Center.x = cuDFNsys::RandomUniform((T)(-model_L * 0.5),
                                                    (T)(model_L * 0.5),
                                                    ((T)rand() / (T)RAND_MAX));
        verts[i].Center.y = cuDFNsys::RandomUniform((T)(-model_L * 0.5),
                                                    (T)(model_L * 0.5),
                                                    ((T)rand() / (T)RAND_MAX));
        verts[i].Center.z = cuDFNsys::RandomUniform((T)(-model_L * 0.5),
                                                    (T)(model_L * 0.5),
                                                    ((T)rand() / (T)RAND_MAX));

        verts[i].NormalVec = cuDFNsys::MakeVector3((T)cuDFNsys::RandomUniform((cuDFNsys::Vector1<T>)-1.0, (cuDFNsys::Vector1<T>)1.0, ((T)rand() / (T)RAND_MAX)),
                                                   (T)cuDFNsys::RandomUniform((cuDFNsys::Vector1<T>)-1.0, (cuDFNsys::Vector1<T>)1.0, ((T)rand() / (T)RAND_MAX)),
                                                   (T)0);
        cuDFNsys::Vector1<T> R_xy = sqrt(verts[i].NormalVec.x * verts[i].NormalVec.x +
                                         verts[i].NormalVec.y * verts[i].NormalVec.y);
        verts[i].NormalVec.z = R_xy / tan(cuDFNsys::RandomFisher(((T)rand() / (T)RAND_MAX),
                                                                 (T)kappa));
        cuDFNsys::Vector1<T> norm_f = sqrt(verts[i].NormalVec.x * verts[i].NormalVec.x +
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
                verts_3D_ptr[(j + 1) % 3] = cuDFNsys::RandomUniform((cuDFNsys::Vector1<T>)-1.0, (cuDFNsys::Vector1<T>)1.0, ((T)rand() / (T)RAND_MAX));
                verts_3D_ptr[(j + 2) % 3] = cuDFNsys::RandomUniform((cuDFNsys::Vector1<T>)-1.0, (cuDFNsys::Vector1<T>)1.0, ((T)rand() / (T)RAND_MAX));
                verts_3D_ptr[j] = -1.0 * (verts_3D_ptr[(j + 1) % 3] * normal_fff[(j + 1) % 3] + verts_3D_ptr[(j + 2) % 3] * normal_fff[(j + 2) % 3]) / normal_fff[j];
                break;
            }
        }

        cuDFNsys::Vector1<T> norm_vert1 = sqrt(verts[i].Verts3D[0].x * verts[i].Verts3D[0].x +
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

        bool gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 0, 1);
        verts[i].ConnectModelSurf[0] = gh;

        gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 0, -1);
        verts[i].ConnectModelSurf[1] = gh;

        gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 1, 1);
        verts[i].ConnectModelSurf[2] = gh;

        gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 1, -1);
        verts[i].ConnectModelSurf[3] = gh;

        gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 2, 1);
        verts[i].ConnectModelSurf[4] = gh;

        gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 2, -1);
        verts[i].ConnectModelSurf[5] = gh;
    }
}; // FracturesCPU
template cuDFNsys::FracturesCPU<double>::FracturesCPU(thrust::host_vector<cuDFNsys::Fracture<double>> &verts,
                                                      unsigned long seed,
                                                      int count,
                                                      double model_L,
                                                      uint ModeSizeDistri,                      // 0 = power law; 1 = lognormal; 2 = uniform; 3 = monosize
                                                      cuDFNsys::Vector4<double> ParaSizeDistri, // when mode = 1, ;
                                                      double kappa,
                                                      double conductivity_powerlaw_exponent, uint Nproc);
template cuDFNsys::FracturesCPU<float>::FracturesCPU(thrust::host_vector<cuDFNsys::Fracture<float>> &verts,
                                                     unsigned long seed,
                                                     int count,
                                                     float model_L,
                                                     uint ModeSizeDistri,                     // 0 = power law; 1 = lognormal; 2 = uniform; 3 = monosize
                                                     cuDFNsys::Vector4<float> ParaSizeDistri, // when mode = 1, ;
                                                     float kappa,
                                                     float conductivity_powerlaw_exponent, uint Nproc);

// ====================================================
// NAME:        FracturesCrossedVertical
// DESCRIPTION: two crossed vertical fractures
// AUTHOR:      Tingchang YIN
// DATE:        20/04/2022
// ====================================================
template <typename T>
__global__ void cuDFNsys::FracturesCrossedVertical(cuDFNsys::Fracture<T> *verts,
                                                   unsigned long seed,
                                                   int count,
                                                   T model_L)
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
        verts[i].Verts3D[0] = cuDFNsys::MakeVector3((T)(-0.5 * model_L),
                                                    (T)0.0,
                                                    (T)(0.5 * model_L));
    else
        verts[i].Verts3D[0] = cuDFNsys::MakeVector3((T)0.0,
                                                    (T)(-0.5 * model_L),
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

    bool gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 0, 1);
    verts[i].ConnectModelSurf[0] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 0, -1);
    verts[i].ConnectModelSurf[1] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 1, 1);
    verts[i].ConnectModelSurf[2] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 1, -1);
    verts[i].ConnectModelSurf[3] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 2, 1);
    verts[i].ConnectModelSurf[4] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 2, -1);
    verts[i].ConnectModelSurf[5] = gh;
}; // FracturesCrossedVertical
template __global__ void cuDFNsys::FracturesCrossedVertical<double>(cuDFNsys::Fracture<double> *verts,
                                                                    unsigned long seed,
                                                                    int count,
                                                                    double model_L);
template __global__ void cuDFNsys::FracturesCrossedVertical<float>(cuDFNsys::Fracture<float> *verts,
                                                                   unsigned long seed,
                                                                   int count,
                                                                   float model_L);

// ====================================================
// NAME:        FracturesBeta50Beta60
// DESCRIPTION: two inclined fractures, with two beta values
// AUTHOR:      Tingchang YIN
// DATE:        20/04/2022
// ====================================================
template <typename T>
__global__ void cuDFNsys::FracturesBeta50Beta60(cuDFNsys::Fracture<T> *verts,
                                                unsigned long seed,
                                                int count,
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
            verts_3D_ptr[(j + 1) % 3] = cuDFNsys::RandomUniform((T)-1.0, (T)1.0, (T)curand_uniform(&state));
            verts_3D_ptr[(j + 2) % 3] = cuDFNsys::RandomUniform((T)-1.0, (T)1.0, (T)curand_uniform(&state));
            verts_3D_ptr[j] = -1.0 * (verts_3D_ptr[(j + 1) % 3] * normal_fff[(j + 1) % 3] + verts_3D_ptr[(j + 2) % 3] * normal_fff[(j + 2) % 3]) / normal_fff[j];
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

    bool gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 0, 1);
    verts[i].ConnectModelSurf[0] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 0, -1);
    verts[i].ConnectModelSurf[1] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 1, 1);
    verts[i].ConnectModelSurf[2] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 1, -1);
    verts[i].ConnectModelSurf[3] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 2, 1);
    verts[i].ConnectModelSurf[4] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 2, -1);
    verts[i].ConnectModelSurf[5] = gh;
}; // FracturesBeta50Beta60
template __global__ void cuDFNsys::FracturesBeta50Beta60<double>(cuDFNsys::Fracture<double> *verts,
                                                                 unsigned long seed,
                                                                 int count,
                                                                 double model_L);
template __global__ void cuDFNsys::FracturesBeta50Beta60<float>(cuDFNsys::Fracture<float> *verts,
                                                                unsigned long seed,
                                                                int count,
                                                                float model_L);

// ====================================================
// NAME:        FracturesIncomplete
// DESCRIPTION: two incomplet fractures
// AUTHOR:      Tingchang YIN
// DATE:        20/04/2022
// ====================================================
template <typename T>
__global__ void cuDFNsys::FracturesIncomplete(cuDFNsys::Fracture<T> *verts,
                                              unsigned long seed,
                                              int count,
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

    bool gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 0, 1);
    verts[i].ConnectModelSurf[0] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 0, -1);
    verts[i].ConnectModelSurf[1] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 1, 1);
    verts[i].ConnectModelSurf[2] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 1, -1);
    verts[i].ConnectModelSurf[3] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 2, 1);
    verts[i].ConnectModelSurf[4] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 2, -1);
    verts[i].ConnectModelSurf[5] = gh;
}; // FracturesIncomplete
template __global__ void cuDFNsys::FracturesIncomplete<double>(cuDFNsys::Fracture<double> *verts,
                                                               unsigned long seed,
                                                               int count,
                                                               double model_L);
template __global__ void cuDFNsys::FracturesIncomplete<float>(cuDFNsys::Fracture<float> *verts,
                                                              unsigned long seed,
                                                              int count,
                                                              float model_L);

// ====================================================
// NAME:        Fractures2DLike
// DESCRIPTION: 2D-like fractures
// AUTHOR:      Tingchang YIN
// DATE:        20/04/2022
// ====================================================
template <typename T>
__global__ void cuDFNsys::Fractures2DLike(cuDFNsys::Fracture<T> *verts,
                                          unsigned long seed,
                                          int count,
                                          T model_L,
                                          T alpha,
                                          T minR,
                                          T maxR)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;

    if (i > count - 1)
        return;

    curandState state;

    curand_init(seed, i, 0, &state);

    T R_ = 0;

    R_ = cuDFNsys::RandomPowerlaw((T)minR, (T)maxR, (T)alpha, (T)curand_uniform(&state));

    verts[i].Radius = R_;
    //printf("%f\n", verts[i].Radius);

    verts[i].NumVertsTruncated = 4;
    //printf("conductivity_powerlaw_exponent: %f, kappa: %f\n", conductivity_powerlaw_exponent, kappa);

    verts[i].Conductivity = (T)(pow(1.0e-3, 3.0) / 12.0);

    verts[i].Center.x = cuDFNsys::RandomUniform((T)(-model_L * 0.5),
                                                (T)(model_L * 0.5),
                                                (T)(curand_uniform(&state)));
    verts[i].Center.y = 0;
    verts[i].Center.z = cuDFNsys::RandomUniform((T)(-model_L * 0.5),
                                                (T)(model_L * 0.5),
                                                (T)(curand_uniform(&state)));

    verts[i].NormalVec = cuDFNsys::MakeVector3((T)cuDFNsys::RandomUniform((T)-1.0, (T)1.0, (T)curand_uniform(&state)),
                                               (T)0,
                                               (T)cuDFNsys::RandomUniform((T)-1.0, (T)1.0, (T)curand_uniform(&state)));
    T norm_f = sqrt(verts[i].NormalVec.x * verts[i].NormalVec.x +
                    verts[i].NormalVec.y * verts[i].NormalVec.y +
                    verts[i].NormalVec.z * verts[i].NormalVec.z);
    verts[i].NormalVec.x /= norm_f;
    verts[i].NormalVec.y /= norm_f;
    verts[i].NormalVec.z /= norm_f;

    //--------------
    cuDFNsys::Vector3<T> rotate_axis = cuDFNsys::MakeVector3(-verts[i].NormalVec.y,
                                                             verts[i].NormalVec.x,
                                                             (T)0);
    T rotate_axis_norm = sqrt(rotate_axis.x * rotate_axis.x +
                              rotate_axis.y * rotate_axis.y +
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

        verts[i].Verts3D[j].y = (verts[i].Verts3D[j].y /
                                 abs(verts[i].Verts3D[j].y)) *
                                0.5 * model_L;
        verts[i].Verts3DTruncated[j].x = verts[i].Verts3D[j].x;
        verts[i].Verts3DTruncated[j].y = verts[i].Verts3D[j].y;
        verts[i].Verts3DTruncated[j].z = verts[i].Verts3D[j].z;
    };

    bool gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 0, 1);
    verts[i].ConnectModelSurf[0] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 0, -1);
    verts[i].ConnectModelSurf[1] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 1, 1);
    verts[i].ConnectModelSurf[2] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 1, -1);
    verts[i].ConnectModelSurf[3] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 2, 1);
    verts[i].ConnectModelSurf[4] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 2, -1);
    verts[i].ConnectModelSurf[5] = gh;
}; // Fractures2DLike
template __global__ void cuDFNsys::Fractures2DLike<double>(cuDFNsys::Fracture<double> *verts,
                                                           unsigned long seed,
                                                           int count,
                                                           double model_L,
                                                           double alpha,
                                                           double minR,
                                                           double maxR);
template __global__ void cuDFNsys::Fractures2DLike<float>(cuDFNsys::Fracture<float> *verts,
                                                          unsigned long seed,
                                                          int count,
                                                          float model_L,
                                                          float alpha,
                                                          float minR,
                                                          float maxR);

// ====================================================
// NAME:        FracturesFour
// DESCRIPTION: Four fractures to verify particle tracking
// AUTHOR:      Tingchang YIN
// DATE:        01/09/2022
// ====================================================
template <typename T>
__global__ void cuDFNsys::FracturesFour(cuDFNsys::Fracture<T> *verts,
                                        unsigned long seed,
                                        int count,
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
        verts[i].Center = cuDFNsys::MakeVector3((T)0., (T)0., (T)0.25 * model_L);
    else if (i == 2)
        verts[i].Center = cuDFNsys::MakeVector3((T)0.25 * model_L, (T)0., (T)(-0.25) * model_L);
    else if (i == 3)
        verts[i].Center = cuDFNsys::MakeVector3((T)(-0.25) * model_L, (T)0., (T)(-0.25) * model_L);

    if (i == 0)
        verts[i].NormalVec = cuDFNsys::MakeVector3((T)0., (T)0., (T)1.); //
    else if (i == 1)
        verts[i].NormalVec = cuDFNsys::MakeVector3(-(T)1., (T)0., (T)0.);
    else if (i == 2)
        verts[i].NormalVec = cuDFNsys::MakeVector3(-(T)1., (T)0., (T)0.);
    else if (i == 3)
        verts[i].NormalVec = cuDFNsys::MakeVector3(-(T)1., (T)0., (T)0.);

    if (i == 0)
        verts[i].Verts3D[0] = cuDFNsys::MakeVector3((T)0.5 * model_L, (T)0.5 * model_L, (T)0.); //
    else
        verts[i].Verts3D[0] = cuDFNsys::MakeVector3((T)0., (T)0.5 * model_L, (T)(0.5) * model_L);

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

    bool gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 0, 1);
    verts[i].ConnectModelSurf[0] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 0, -1);
    verts[i].ConnectModelSurf[1] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 1, 1);
    verts[i].ConnectModelSurf[2] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 1, -1);
    verts[i].ConnectModelSurf[3] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 2, 1);
    verts[i].ConnectModelSurf[4] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 2, -1);
    verts[i].ConnectModelSurf[5] = gh;
};
template __global__ void cuDFNsys::FracturesFour<double>(cuDFNsys::Fracture<double> *verts,
                                                         unsigned long seed,
                                                         int count,
                                                         double model_L);
template __global__ void cuDFNsys::FracturesFour<float>(cuDFNsys::Fracture<float> *verts,
                                                        unsigned long seed,
                                                        int count,
                                                        float model_L);

// ====================================================
// NAME:        FracturesChangeDomainSize
// DESCRIPTION: Change Domain size
// AUTHOR:      Tingchang YIN
// DATE:        13/09/2022
// ====================================================
template <typename T>
__global__ void cuDFNsys::FracturesChangeDomainSize(cuDFNsys::Fracture<T> *verts,
                                                    int count,
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

    bool gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 0, 1);
    verts[i].ConnectModelSurf[0] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 0, -1);
    verts[i].ConnectModelSurf[1] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 1, 1);
    verts[i].ConnectModelSurf[2] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 1, -1);
    verts[i].ConnectModelSurf[3] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 2, 1);
    verts[i].ConnectModelSurf[4] = gh;

    gh = cuDFNsys::TruncateFracture<T>(&verts[i], model_L, 2, -1);
    verts[i].ConnectModelSurf[5] = gh;
}; // FracturesChangeDomainSize
template __global__ void cuDFNsys::FracturesChangeDomainSize<double>(cuDFNsys::Fracture<double> *verts,
                                                                     int count,
                                                                     double model_L);
template __global__ void cuDFNsys::FracturesChangeDomainSize<float>(cuDFNsys::Fracture<float> *verts,
                                                                    int count,
                                                                    float model_L);