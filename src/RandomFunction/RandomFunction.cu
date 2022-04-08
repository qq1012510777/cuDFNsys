#include "RandomFunction/RandomFunction.cuh"

// ====================================================
// NAME:        RandomPowerlaw
// DESCRIPTION: return a value having a power law distribution.
// AUTHOR:      Tingchang YIN
// DATE:        04/04/2022
// ====================================================
__device__ __host__ float cuDFNsys::RandomPowerlaw(float x0,
                                                   float x1,
                                                   float alpha_g,
                                                   float rand_0_1)
{
    float x_g = (pow(x1, 1.0f - alpha_g) - pow(x0, 1.0f - alpha_g)) * rand_0_1 + pow(x0, 1.0f - alpha_g);
    x_g = pow(x_g, 1.0f / (1.0f - alpha_g));
    return x_g;
}; // RandomPowerlaw

// ====================================================
// NAME:        RandomUniform
// DESCRIPTION: return a value having a uniform distribution.
// AUTHOR:      Tingchang YIN
// DATE:        04/04/2022
// ====================================================
__device__ __host__ float cuDFNsys::RandomUniform(float a,
                                                  float b,
                                                  float rand_0_1)
{
    float random = a + (b - a) * rand_0_1;
    return random;
}; // RandomUniform

// ====================================================
// NAME:        RandomFisher
// DESCRIPTION: return a radian value having a fisher
//              distribution.
// AUTHOR:      Tingchang YIN
// DATE:        04/04/2022
// ====================================================
__device__ __host__ float cuDFNsys::RandomFisher(float rand_0_1,
                                                 float fisher_k)
{
    if (fisher_k == 0.0f)
        return (rand_0_1 * 0.5f * M_PI);

UKpo100:;
    float theta = 0.0f;

    if (fisher_k < 88.0f)
        theta = acos(log(exp(fisher_k) - (exp(fisher_k) - exp(-1.0 * fisher_k)) * rand_0_1) / fisher_k);
    else
    {
        double kappa = (double)fisher_k;
        double thetay = acos(log(exp(kappa) - (exp(kappa) - exp(-1.0 * kappa)) * (double)rand_0_1) / kappa);
        theta = (float)thetay;
    }

    if (theta > 0.5f * M_PI)
        theta = M_PI - theta;

    if (theta == 0.5f * M_PI)
        goto UKpo100;

    return theta;
}; // RandomFisher