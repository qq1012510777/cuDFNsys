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

#include "RandomFunction/RandomFunction.cuh"

// ====================================================
// NAME:        RandomPowerlaw
// DESCRIPTION: return a value having a power law distribution.
// AUTHOR:      Tingchang YIN
// DATE:        04/04/2022
// ====================================================
template <typename T>
__device__ __host__ T cuDFNsys::RandomPowerlaw(T x0,
                                               T x1,
                                               T alpha_g,
                                               T rand_0_1)
{
    T x_g = (pow(x1, 1.0f - alpha_g) - pow(x0, 1.0f - alpha_g)) * rand_0_1 + pow(x0, 1.0f - alpha_g);
    x_g = pow(x_g, 1.0f / (1.0f - alpha_g));
    return x_g;
}; // RandomPowerlaw
template __device__ __host__ double cuDFNsys::RandomPowerlaw(double x0,
                                                             double x1,
                                                             double alpha_g,
                                                             double rand_0_1);
template __device__ __host__ float cuDFNsys::RandomPowerlaw(float x0,
                                                            float x1,
                                                            float alpha_g,
                                                            float rand_0_1);

// ====================================================
// NAME:        RandomUniform
// DESCRIPTION: return a value having a uniform distribution.
// AUTHOR:      Tingchang YIN
// DATE:        04/04/2022
// ====================================================
template <typename T>
__device__ __host__ T cuDFNsys::RandomUniform(T a,
                                              T b,
                                              T rand_0_1)
{
    T random = a + (b - a) * rand_0_1;
    return random;
}; // RandomUniform
template __device__ __host__ double cuDFNsys::RandomUniform(double a,
                                                            double b,
                                                            double rand_0_1);
template __device__ __host__ float cuDFNsys::RandomUniform(float a,
                                                           float b,
                                                           float rand_0_1);

// ====================================================
// NAME:        RandomFisher
// DESCRIPTION: return a radian value having a fisher
//              distribution.
// AUTHOR:      Tingchang YIN
// DATE:        04/04/2022
// ====================================================
template <typename T>
__device__ __host__ T cuDFNsys::RandomFisher(T rand_0_1,
                                             T fisher_k)
{
    if (fisher_k == (T)0.0)
        return (rand_0_1 * (T)0.5 * M_PI);

UKpo100:;
    T theta = 0.0;

    if (fisher_k < (T)88.0)
        theta = acos(log(exp(fisher_k) - (exp(fisher_k) - exp(-1.0 * fisher_k)) * rand_0_1) / fisher_k);
    else
    {
        double kappa = (double)fisher_k;
        double thetay = acos(log(exp(kappa) - (exp(kappa) - exp(-1.0 * kappa)) * (double)rand_0_1) / kappa);
        theta = (T)thetay;
    }

    if (theta > (T)(0.5 * M_PI))
        theta = M_PI - theta;

    if (theta == (T)(0.5 * M_PI))
        goto UKpo100;

    return theta;
}; // RandomFisher
template __device__ __host__ double cuDFNsys::RandomFisher(double rand_0_1,
                                                           double fisher_k);
template __device__ __host__ float cuDFNsys::RandomFisher(float rand_0_1,
                                                          float fisher_k);
// ====================================================
// NAME:        RandomLognormal
// DESCRIPTION: return a value having lognormal
//              distribution
// AUTHOR:      Tingchang YIN
// DATE:        14/04/2022
// ====================================================
template <typename T>
__device__ __host__ T cuDFNsys::RandomLognormal(T mean,
                                                T variance,
                                                T min,
                                                T max,
                                                T rand_0_1)
{
    T mu_1 = log(mean * mean / (pow(variance + mean * mean, 0.5f)));
    T sigma_1 = pow(log(1 + ((T)variance) / (mean * mean)), 0.5f);
    //sigma_1 is input std. deviation

    T l = (log(min) - mu_1) / sigma_1;
    T u = (log(max) - mu_1) / sigma_1;

    T p_l = cuDFNsys::StandardNormalCDF(l);
    T p_u = cuDFNsys::StandardNormalCDF(u);

    T x_prime = p_l + (p_u - p_l) * rand_0_1;

    T z = cuDFNsys::StandardNormalCDFInv(x_prime);

    T value = exp(z * sigma_1 + mu_1);

    return value;
};
template __device__ __host__ double cuDFNsys::RandomLognormal(double mean,
                                                              double variance,
                                                              double min,
                                                              double max,
                                                              double rand_0_1);
template __device__ __host__ float cuDFNsys::RandomLognormal(float mean,
                                                             float variance,
                                                             float min,
                                                             float max,
                                                             float rand_0_1);

// ====================================================
// NAME:        StandardNormalCDF
// DESCRIPTION: return CDF of standard normal distribution (-inf, x)
// AUTHOR:      Tingchang YIN
// DATE:        14/04/2022
// ====================================================
template <typename T>
__device__ __host__ T cuDFNsys::StandardNormalCDF(T x)
{
    T a1 = 0.398942280444;
    T a2 = 0.399903438504;
    T a3 = 5.75885480458;
    T a4 = 29.8213557808;
    T a5 = 2.62433121679;
    T a6 = 48.6959930692;
    T a7 = 5.92885724438;
    T b0 = 0.398942280385;
    T b1 = 3.8052e-08;
    T b2 = 1.00000615302;
    T b3 = 3.98064794e-04;
    T b4 = 1.98615381364;
    T b5 = 0.151679116635;
    T b6 = 5.29330324926;
    T b7 = 4.8385912808;
    T b8 = 15.1508972451;
    T b9 = 0.742380924027;
    T b10 = 30.789933034;
    T b11 = 3.99019417011;
    T cdf;
    T q;
    T y;
    //
    //  |X| <= 1.28.
    //
    if (fabs(x) <= (T)1.28)
    {
        y = (T)0.5 * x * x;

        q = (T)0.5 - fabs(x) * (a1 - a2 * y / (y + a3 - a4 / (y + a5 + a6 / (y + a7))));
        //
        //  1.28 < |X| <= 12.7
        //
    }
    else if (fabs(x) <= (T)12.7)
    {
        y = (T)0.5 * x * x;

        q = exp(-y) * b0 / (fabs(x) - b1 + b2 / (fabs(x) + b3 + b4 / (fabs(x) - b5 + b6 / (fabs(x) + b7 - b8 / (fabs(x) + b9 + b10 / (fabs(x) + b11))))));
        //
        //  12.7 < |X|
        //
    }
    else
    {
        q = (T)0.0;
    }
    //
    //  Take account of negative X.
    //
    if (x < (T)0.0)
    {
        cdf = q;
    }
    else
    {
        cdf = (T)1.0 - q;
    }

    return cdf;
};
template __device__ __host__ double cuDFNsys::StandardNormalCDF(double x);
template __device__ __host__ float cuDFNsys::StandardNormalCDF(float x);

// ====================================================
// NAME:        StandardNormalCDFInv
// DESCRIPTION: return a value that is calculated by
//              inverting the standard normal CDF.
// AUTHOR:      Tingchang YIN
// DATE:        15/04/2022
// ====================================================
template <typename T>
__device__ __host__ T cuDFNsys::StandardNormalCDFInv(T p)
{
    T a[8] = {
        3.3871328727963666080, 1.3314166789178437745e+2,
        1.9715909503065514427e+3, 1.3731693765509461125e+4,
        4.5921953931549871457e+4, 6.7265770927008700853e+4,
        3.3430575583588128105e+4, 2.5090809287301226727e+3};
    T b[8] = {
        1.0, 4.2313330701600911252e+1,
        6.8718700749205790830e+2, 5.3941960214247511077e+3,
        2.1213794301586595867e+4, 3.9307895800092710610e+4,
        2.8729085735721942674e+4, 5.2264952788528545610e+3};
    T c[8] = {
        1.42343711074968357734, 4.63033784615654529590,
        5.76949722146069140550, 3.64784832476320460504,
        1.27045825245236838258, 2.41780725177450611770e-1,
        2.27238449892691845833e-2, 7.74545014278341407640e-4};
    T const1 = 0.180625;
    T const2 = 1.6;
    T d[8] = {
        1.0, 2.05319162663775882187,
        1.67638483018380384940, 6.89767334985100004550e-1,
        1.48103976427480074590e-1, 1.51986665636164571966e-2,
        5.47593808499534494600e-4, 1.05075007164441684324e-9};
    T e[8] = {
        6.65790464350110377720, 5.46378491116411436990,
        1.78482653991729133580, 2.96560571828504891230e-1,
        2.65321895265761230930e-2, 1.24266094738807843860e-3,
        2.71155556874348757815e-5, 2.01033439929228813265e-7};
    T f[8] = {
        1.0, 5.99832206555887937690e-1,
        1.36929880922735805310e-1, 1.48753612908506148525e-2,
        7.86869131145613259100e-4, 1.84631831751005468180e-5,
        1.42151175831644588870e-7, 2.04426310338993978564e-15};
    T q;
    T r;
    const T r8_huge = 1.0E+30;
    T split1 = 0.425;
    T split2 = 5.0;
    T value;

    if (p <= 0.0f)
    {
        value = -r8_huge;
        return value;
    }

    if (1.0f <= p)
    {
        value = r8_huge;
        return value;
    }

    q = p - 0.5f;

    if (fabs(q) <= split1)
    {
        r = const1 - q * q;
        value = q * cuDFNsys::R8PolyValueHorner(7, a, r) / cuDFNsys::R8PolyValueHorner(7, b, r);
    }
    else
    {
        if (q < 0.0f)
        {
            r = p;
        }
        else
        {
            r = 1.0f - p;
        }

        if (r <= 0.0f)
        {
            value = r8_huge;
        }
        else
        {
            r = sqrt(-log(r));

            if (r <= split2)
            {
                r = r - const2;
                value = cuDFNsys::R8PolyValueHorner(7, c, r) / cuDFNsys::R8PolyValueHorner(7, d, r);
            }
            else
            {
                r = r - split2;
                value = cuDFNsys::R8PolyValueHorner(7, e, r) / cuDFNsys::R8PolyValueHorner(7, f, r);
            }
        }

        if (q < 0.0f)
        {
            value = -value;
        }
    }

    return value;
};
template __device__ __host__ double cuDFNsys::StandardNormalCDFInv(double p);
template __device__ __host__ float cuDFNsys::StandardNormalCDFInv(float p);

// ====================================================
// NAME:        R8PolyValueHorner
// DESCRIPTION: A stable algorithm to calculate
//              polynomial with 8th degrees
// AUTHOR:      Tingchang YIN
// DATE:        15/04/2022
// ====================================================
template <typename T>
__device__ __host__ T cuDFNsys::R8PolyValueHorner(int m, T c[], T x)
{
    int i;
    T value;

    value = c[m];

    for (i = m - 1; 0 <= i; i--)
    {
        value = value * x + c[i];
    }

    return value;
};
template __device__ __host__ double cuDFNsys::R8PolyValueHorner(int m, double c[], double x);
template __device__ __host__ float cuDFNsys::R8PolyValueHorner(int m, float c[], float x);

// ====================================================
// NAME:        RandomStandardNormal
// DESCRIPTION: from variable following uniform distribution of (0, 1)
//              to value following standard normal distribution
// AUTHOR:      Tingchang YIN
// DATE:        10/11/2022
// ====================================================
template <typename T>
__device__ __host__ T cuDFNsys::RandomStandardNormal(T rand_0_1_1, T rand_0_1_2)
{
    return sqrt(-2.0 * log(rand_0_1_1)) * cos(2.0 * M_PI * rand_0_1_2);
}; // RandomStandardNormal
template __device__ __host__ double cuDFNsys::RandomStandardNormal(double rand_0_1_1, double rand_0_1_2);
template __device__ __host__ float cuDFNsys::RandomStandardNormal(float rand_0_1_1, float rand_0_1_2);