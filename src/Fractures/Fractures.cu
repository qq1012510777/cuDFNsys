#include "Fractures/Fractures.cuh"

// ====================================================
// NAME:        Fractures
// DESCRIPTION: Fractures in a DFN
// AUTHOR:      Tingchang YIN
// DATE:        04/04/2022
// ====================================================
__global__ void cuDFNsys::Fractures(cuDFNsys::Fracture *verts,
                                    unsigned long seed,
                                    int count,
                                    float model_L,
                                    uint ModeSizeDistri,
                                    float4 ParaSizeDistri,
                                    float kappa,
                                    float conductivity_powerlaw_exponent)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;

    if (i > count - 1)
        return;

    curandState state;

    curand_init(seed, i, 0, &state);

    float R_ = 0;

    // if (alpha == 0 && abs(minR - maxR) < 1e-7)
    //     R_ = minR;
    // else if (alpha == 0 && abs(minR - maxR) > 1e-7)
    //     R_ = cuDFNsys::RandomUniform(minR, maxR, curand_uniform(&state));
    // else
    //     R_ = cuDFNsys::RandomPowerlaw(minR, maxR, alpha, curand_uniform(&state));

    if (ModeSizeDistri == 0)
        R_ = cuDFNsys::RandomPowerlaw(ParaSizeDistri.y, ParaSizeDistri.z,
                                      ParaSizeDistri.x, curand_uniform(&state));
    else if (ModeSizeDistri == 1)
        R_ = cuDFNsys::RandomLognormal(ParaSizeDistri.x,
                                       ParaSizeDistri.y,
                                       ParaSizeDistri.z,
                                       ParaSizeDistri.w, curand_uniform(&state));
    else if (ModeSizeDistri == 2)
        R_ = cuDFNsys::RandomUniform(ParaSizeDistri.x,
                                     ParaSizeDistri.y, curand_uniform(&state));
    else if (ModeSizeDistri == 3)
        R_ = ParaSizeDistri.x;

    verts[i].Radius = R_;
    //printf("%f\n", verts[i].Radius);

    verts[i].NumVertsTruncated = 4;
    //printf("conductivity_powerlaw_exponent: %f, kappa: %f\n", conductivity_powerlaw_exponent, kappa);

    if (conductivity_powerlaw_exponent == 0)
        verts[i].Conductivity = 1.0;
    else
        verts[i].Conductivity = (1.0e-11) * pow(R_, 3.0 * conductivity_powerlaw_exponent);

    verts[i].Center.x = cuDFNsys::RandomUniform(-model_L * 0.5, model_L * 0.5,
                                                curand_uniform(&state));
    verts[i].Center.y = cuDFNsys::RandomUniform(-model_L * 0.5, model_L * 0.5,
                                                curand_uniform(&state));
    verts[i].Center.z = cuDFNsys::RandomUniform(-model_L * 0.5, model_L * 0.5,
                                                curand_uniform(&state));

    verts[i].NormalVec = make_float3(cuDFNsys::RandomUniform(-1.0, 1.0, curand_uniform(&state)),
                                     cuDFNsys::RandomUniform(-1.0, 1.0, curand_uniform(&state)),
                                     0);
    float R_xy = sqrt(verts[i].NormalVec.x * verts[i].NormalVec.x +
                      verts[i].NormalVec.y * verts[i].NormalVec.y);
    verts[i].NormalVec.z = R_xy / tan(cuDFNsys::RandomFisher(curand_uniform(&state), kappa));
    float norm_f = sqrt(verts[i].NormalVec.x * verts[i].NormalVec.x +
                        verts[i].NormalVec.y * verts[i].NormalVec.y +
                        verts[i].NormalVec.z * verts[i].NormalVec.z);
    verts[i].NormalVec.x /= norm_f;
    verts[i].NormalVec.y /= norm_f;
    verts[i].NormalVec.z /= norm_f;

    float *normal_fff = &verts[i].NormalVec.x;
    float *verts_3D_ptr = &verts[i].Verts3D[0].x;
    for (int j = 0; j < 3; ++j)
    {
        if (abs(normal_fff[j]) > 1e-3)
        {
            verts_3D_ptr[(j + 1) % 3] = cuDFNsys::RandomUniform(-1.0, 1.0, curand_uniform(&state));
            verts_3D_ptr[(j + 2) % 3] = cuDFNsys::RandomUniform(-1.0, 1.0, curand_uniform(&state));
            verts_3D_ptr[j] = -1.0 * (verts_3D_ptr[(j + 1) % 3] * normal_fff[(j + 1) % 3] + verts_3D_ptr[(j + 2) % 3] * normal_fff[(j + 2) % 3]) / normal_fff[j];
            break;
        }
    }

    float norm_vert1 = sqrt(verts[i].Verts3D[0].x * verts[i].Verts3D[0].x +
                            verts[i].Verts3D[0].y * verts[i].Verts3D[0].y +
                            verts[i].Verts3D[0].z * verts[i].Verts3D[0].z);
    norm_vert1 = R_ / norm_vert1;
    verts[i].Verts3D[0].x *= norm_vert1;
    verts[i].Verts3D[0].y *= norm_vert1;
    verts[i].Verts3D[0].z *= norm_vert1;
    verts[i].Verts3D[2].x = -1.0 * verts[i].Verts3D[0].x;
    verts[i].Verts3D[2].y = -1.0 * verts[i].Verts3D[0].y;
    verts[i].Verts3D[2].z = -1.0 * verts[i].Verts3D[0].z;

    verts[i].Verts3D[1] = cuDFNsys::CrossProductFloat3(verts[i].NormalVec,
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

    bool gh = cuDFNsys::TruncateFracture(&verts[i], model_L, 0, 1);
    verts[i].ConnectModelSurf[0] = gh;

    gh = cuDFNsys::TruncateFracture(&verts[i], model_L, 0, -1);
    verts[i].ConnectModelSurf[1] = gh;

    gh = cuDFNsys::TruncateFracture(&verts[i], model_L, 1, 1);
    verts[i].ConnectModelSurf[2] = gh;

    gh = cuDFNsys::TruncateFracture(&verts[i], model_L, 1, -1);
    verts[i].ConnectModelSurf[3] = gh;

    gh = cuDFNsys::TruncateFracture(&verts[i], model_L, 2, 1);
    verts[i].ConnectModelSurf[4] = gh;

    gh = cuDFNsys::TruncateFracture(&verts[i], model_L, 2, -1);
    verts[i].ConnectModelSurf[5] = gh;
}; // Fractures

// ====================================================
// NAME:        FracturesCrossedVertical
// DESCRIPTION: two crossed vertical fractures
// AUTHOR:      Tingchang YIN
// DATE:        20/04/2022
// ====================================================
__global__ void cuDFNsys::FracturesCrossedVertical(cuDFNsys::Fracture *verts,
                                                   unsigned long seed,
                                                   int count,
                                                   float model_L)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;

    if (i >= 2)
        return;

    curandState state;

    curand_init(seed, i, 0, &state);

    float R_ = sqrt(model_L * model_L / 2.0f);

    verts[i].Radius = R_;
    //printf("%f\n", verts[i].Radius);

    verts[i].NumVertsTruncated = 4;
    //printf("conductivity_powerlaw_exponent: %f, kappa: %f\n", conductivity_powerlaw_exponent, kappa);

    verts[i].Conductivity = 1;

    verts[i].Center.x = 0;
    verts[i].Center.y = 0;
    verts[i].Center.z = 0;

    if (i == 0)
        verts[i].NormalVec = make_float3(0, 1.0f, 0);
    else
        verts[i].NormalVec = make_float3(1.0f, 0, 0);

    if (i == 0)
        verts[i].Verts3D[0] = make_float3(-0.5 * model_L, 0, 0.5 * model_L);
    else
        verts[i].Verts3D[0] = make_float3(0, -0.5 * model_L, 0.5 * model_L);

    float norm_vert1 = sqrt(verts[i].Verts3D[0].x * verts[i].Verts3D[0].x +
                            verts[i].Verts3D[0].y * verts[i].Verts3D[0].y +
                            verts[i].Verts3D[0].z * verts[i].Verts3D[0].z);
    norm_vert1 = R_ / norm_vert1;
    verts[i].Verts3D[0].x *= norm_vert1;
    verts[i].Verts3D[0].y *= norm_vert1;
    verts[i].Verts3D[0].z *= norm_vert1;
    verts[i].Verts3D[2].x = -1.0 * verts[i].Verts3D[0].x;
    verts[i].Verts3D[2].y = -1.0 * verts[i].Verts3D[0].y;
    verts[i].Verts3D[2].z = -1.0 * verts[i].Verts3D[0].z;

    verts[i].Verts3D[1] = cuDFNsys::CrossProductFloat3(verts[i].NormalVec,
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

    bool gh = cuDFNsys::TruncateFracture(&verts[i], model_L, 0, 1);
    verts[i].ConnectModelSurf[0] = gh;

    gh = cuDFNsys::TruncateFracture(&verts[i], model_L, 0, -1);
    verts[i].ConnectModelSurf[1] = gh;

    gh = cuDFNsys::TruncateFracture(&verts[i], model_L, 1, 1);
    verts[i].ConnectModelSurf[2] = gh;

    gh = cuDFNsys::TruncateFracture(&verts[i], model_L, 1, -1);
    verts[i].ConnectModelSurf[3] = gh;

    gh = cuDFNsys::TruncateFracture(&verts[i], model_L, 2, 1);
    verts[i].ConnectModelSurf[4] = gh;

    gh = cuDFNsys::TruncateFracture(&verts[i], model_L, 2, -1);
    verts[i].ConnectModelSurf[5] = gh;
}; // FracturesCrossedVertical

// ====================================================
// NAME:        FracturesBeta50Beta60
// DESCRIPTION: two inclined fractures, with two beta values
// AUTHOR:      Tingchang YIN
// DATE:        20/04/2022
// ====================================================
__global__ void cuDFNsys::FracturesBeta50Beta60(cuDFNsys::Fracture *verts,
                                                unsigned long seed,
                                                int count,
                                                float model_L)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;

    if (i >= 2)
        return;

    curandState state;

    curand_init(seed, i, 0, &state);

    float R_ = model_L * 3.0f;

    verts[i].Radius = R_;
    //printf("%f\n", verts[i].Radius);

    verts[i].NumVertsTruncated = 4;
    //printf("conductivity_powerlaw_exponent: %f, kappa: %f\n", conductivity_powerlaw_exponent, kappa);

    verts[i].Conductivity = 1.0;

    if (i == 0)
        verts[i].NormalVec = make_float3(1, 0, 0);
    else if (i == 1)
        verts[i].NormalVec = make_float3(-1, 0, 0);

    float R_xy = sqrt(verts[i].NormalVec.x * verts[i].NormalVec.x +
                      verts[i].NormalVec.y * verts[i].NormalVec.y);
    if (i == 0)
        verts[i].NormalVec.z = R_xy / tan(50.0 / 180.0 * M_PI);
    else if (i == 1)
        verts[i].NormalVec.z = R_xy / tan(1.0 / 3.0 * M_PI);
    float norm_f = sqrt(verts[i].NormalVec.x * verts[i].NormalVec.x +
                        verts[i].NormalVec.y * verts[i].NormalVec.y +
                        verts[i].NormalVec.z * verts[i].NormalVec.z);
    verts[i].NormalVec.x /= norm_f;
    verts[i].NormalVec.y /= norm_f;
    verts[i].NormalVec.z /= norm_f;

    verts[i].Center = make_float3(0, 0, 0);

    //----------
    float *normal_fff = &verts[i].NormalVec.x;
    float *verts_3D_ptr = &verts[i].Verts3D[0].x;
    for (int j = 0; j < 3; ++j)
    {
        if (abs(normal_fff[j]) > 1e-3)
        {
            verts_3D_ptr[(j + 1) % 3] = cuDFNsys::RandomUniform(-1.0, 1.0, curand_uniform(&state));
            verts_3D_ptr[(j + 2) % 3] = cuDFNsys::RandomUniform(-1.0, 1.0, curand_uniform(&state));
            verts_3D_ptr[j] = -1.0 * (verts_3D_ptr[(j + 1) % 3] * normal_fff[(j + 1) % 3] + verts_3D_ptr[(j + 2) % 3] * normal_fff[(j + 2) % 3]) / normal_fff[j];
            break;
        }
    }

    float norm_vert1 = sqrt(verts[i].Verts3D[0].x * verts[i].Verts3D[0].x +
                            verts[i].Verts3D[0].y * verts[i].Verts3D[0].y +
                            verts[i].Verts3D[0].z * verts[i].Verts3D[0].z);
    norm_vert1 = R_ / norm_vert1;
    verts[i].Verts3D[0].x *= norm_vert1;
    verts[i].Verts3D[0].y *= norm_vert1;
    verts[i].Verts3D[0].z *= norm_vert1;
    verts[i].Verts3D[2].x = -1.0 * verts[i].Verts3D[0].x;
    verts[i].Verts3D[2].y = -1.0 * verts[i].Verts3D[0].y;
    verts[i].Verts3D[2].z = -1.0 * verts[i].Verts3D[0].z;

    verts[i].Verts3D[1] = cuDFNsys::CrossProductFloat3(verts[i].NormalVec,
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

    bool gh = cuDFNsys::TruncateFracture(&verts[i], model_L, 0, 1);
    verts[i].ConnectModelSurf[0] = gh;

    gh = cuDFNsys::TruncateFracture(&verts[i], model_L, 0, -1);
    verts[i].ConnectModelSurf[1] = gh;

    gh = cuDFNsys::TruncateFracture(&verts[i], model_L, 1, 1);
    verts[i].ConnectModelSurf[2] = gh;

    gh = cuDFNsys::TruncateFracture(&verts[i], model_L, 1, -1);
    verts[i].ConnectModelSurf[3] = gh;

    gh = cuDFNsys::TruncateFracture(&verts[i], model_L, 2, 1);
    verts[i].ConnectModelSurf[4] = gh;

    gh = cuDFNsys::TruncateFracture(&verts[i], model_L, 2, -1);
    verts[i].ConnectModelSurf[5] = gh;
}; // FracturesBeta50Beta60

// ====================================================
// NAME:        FracturesIncomplete
// DESCRIPTION: two incomplet fractures
// AUTHOR:      Tingchang YIN
// DATE:        20/04/2022
// ====================================================
__global__ void cuDFNsys::FracturesIncomplete(cuDFNsys::Fracture *verts,
                                              unsigned long seed,
                                              int count,
                                              float model_L)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;

    if (i >= 2)
        return;

    curandState state;

    curand_init(seed, i, 0, &state);

    float R_ = (i == 0 ? 0.55 : 0.60) * model_L;
    ;

    verts[i].Radius = R_;
    //printf("%f\n", verts[i].Radius);

    verts[i].NumVertsTruncated = 4;
    //printf("conductivity_powerlaw_exponent: %f, kappa: %f\n", conductivity_powerlaw_exponent, kappa);

    verts[i].Conductivity = 1.0;

    if (i == 0)
        verts[i].NormalVec = make_float3(1, 0, 0);
    else if (i == 1)
        verts[i].NormalVec = make_float3(0, 1, 0);

    verts[i].Center = make_float3(0, 0, 0);

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

    float norm_vert1 = sqrt(verts[i].Verts3D[0].x * verts[i].Verts3D[0].x +
                            verts[i].Verts3D[0].y * verts[i].Verts3D[0].y +
                            verts[i].Verts3D[0].z * verts[i].Verts3D[0].z);
    norm_vert1 = R_ / norm_vert1;
    verts[i].Verts3D[0].x *= norm_vert1;
    verts[i].Verts3D[0].y *= norm_vert1;
    verts[i].Verts3D[0].z *= norm_vert1;
    verts[i].Verts3D[2].x = -1.0 * verts[i].Verts3D[0].x;
    verts[i].Verts3D[2].y = -1.0 * verts[i].Verts3D[0].y;
    verts[i].Verts3D[2].z = -1.0 * verts[i].Verts3D[0].z;

    verts[i].Verts3D[1] = cuDFNsys::CrossProductFloat3(verts[i].NormalVec,
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

    bool gh = cuDFNsys::TruncateFracture(&verts[i], model_L, 0, 1);
    verts[i].ConnectModelSurf[0] = gh;

    gh = cuDFNsys::TruncateFracture(&verts[i], model_L, 0, -1);
    verts[i].ConnectModelSurf[1] = gh;

    gh = cuDFNsys::TruncateFracture(&verts[i], model_L, 1, 1);
    verts[i].ConnectModelSurf[2] = gh;

    gh = cuDFNsys::TruncateFracture(&verts[i], model_L, 1, -1);
    verts[i].ConnectModelSurf[3] = gh;

    gh = cuDFNsys::TruncateFracture(&verts[i], model_L, 2, 1);
    verts[i].ConnectModelSurf[4] = gh;

    gh = cuDFNsys::TruncateFracture(&verts[i], model_L, 2, -1);
    verts[i].ConnectModelSurf[5] = gh;
}; // FracturesIncomplete

// ====================================================
// NAME:        Fractures2DLike
// DESCRIPTION: 2D-like fractures
// AUTHOR:      Tingchang YIN
// DATE:        20/04/2022
// ====================================================
__global__ void cuDFNsys::Fractures2DLike(cuDFNsys::Fracture *verts,
                                          unsigned long seed,
                                          int count,
                                          float model_L,
                                          float alpha,
                                          float minR,
                                          float maxR)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;

    if (i > count - 1)
        return;

    curandState state;

    curand_init(seed, i, 0, &state);

    float R_ = 0;

    R_ = cuDFNsys::RandomPowerlaw(minR, maxR, alpha, curand_uniform(&state));

    verts[i].Radius = R_;
    //printf("%f\n", verts[i].Radius);

    verts[i].NumVertsTruncated = 4;
    //printf("conductivity_powerlaw_exponent: %f, kappa: %f\n", conductivity_powerlaw_exponent, kappa);

    verts[i].Conductivity = 1.0;

    verts[i].Center.x = cuDFNsys::RandomUniform(-model_L * 0.5, model_L * 0.5,
                                                curand_uniform(&state));
    verts[i].Center.y = 0;
    verts[i].Center.z = cuDFNsys::RandomUniform(-model_L * 0.5, model_L * 0.5,
                                                curand_uniform(&state));

    verts[i].NormalVec = make_float3(cuDFNsys::RandomUniform(-1.0, 1.0, curand_uniform(&state)),
                                     0,
                                     cuDFNsys::RandomUniform(-1.0, 1.0, curand_uniform(&state)));
    float norm_f = sqrt(verts[i].NormalVec.x * verts[i].NormalVec.x +
                        verts[i].NormalVec.y * verts[i].NormalVec.y +
                        verts[i].NormalVec.z * verts[i].NormalVec.z);
    verts[i].NormalVec.x /= norm_f;
    verts[i].NormalVec.y /= norm_f;
    verts[i].NormalVec.z /= norm_f;

    //--------------
    float3 rotate_axis = make_float3(-verts[i].NormalVec.y,
                                     verts[i].NormalVec.x,
                                     0);
    float rotate_axis_norm = sqrt(rotate_axis.x * rotate_axis.x +
                                  rotate_axis.y * rotate_axis.y +
                                  rotate_axis.z * rotate_axis.z);
    rotate_axis.x /= rotate_axis_norm;
    rotate_axis.y /= rotate_axis_norm;
    rotate_axis.z /= rotate_axis_norm;

    cuDFNsys::Quaternion qua;
    qua = qua.DescribeRotation(verts[i].NormalVec, 45.0 / 180.0 * M_PI);
    verts[i].Verts3D[0] = qua.Rotate(rotate_axis);

    //--------------
    float norm_vert1 = sqrt(verts[i].Verts3D[0].x * verts[i].Verts3D[0].x +
                            verts[i].Verts3D[0].y * verts[i].Verts3D[0].y +
                            verts[i].Verts3D[0].z * verts[i].Verts3D[0].z);
    norm_vert1 = R_ / norm_vert1;
    verts[i].Verts3D[0].x *= norm_vert1;
    verts[i].Verts3D[0].y *= norm_vert1;
    verts[i].Verts3D[0].z *= norm_vert1;
    verts[i].Verts3D[2].x = -1.0 * verts[i].Verts3D[0].x;
    verts[i].Verts3D[2].y = -1.0 * verts[i].Verts3D[0].y;
    verts[i].Verts3D[2].z = -1.0 * verts[i].Verts3D[0].z;

    verts[i].Verts3D[1] = cuDFNsys::CrossProductFloat3(verts[i].NormalVec,
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

    bool gh = cuDFNsys::TruncateFracture(&verts[i], model_L, 0, 1);
    verts[i].ConnectModelSurf[0] = gh;

    gh = cuDFNsys::TruncateFracture(&verts[i], model_L, 0, -1);
    verts[i].ConnectModelSurf[1] = gh;

    gh = cuDFNsys::TruncateFracture(&verts[i], model_L, 1, 1);
    verts[i].ConnectModelSurf[2] = gh;

    gh = cuDFNsys::TruncateFracture(&verts[i], model_L, 1, -1);
    verts[i].ConnectModelSurf[3] = gh;

    gh = cuDFNsys::TruncateFracture(&verts[i], model_L, 2, 1);
    verts[i].ConnectModelSurf[4] = gh;

    gh = cuDFNsys::TruncateFracture(&verts[i], model_L, 2, -1);
    verts[i].ConnectModelSurf[5] = gh;
}; // Fractures2DLike