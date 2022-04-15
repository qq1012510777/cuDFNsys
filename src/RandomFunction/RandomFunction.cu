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

// ====================================================
// NAME:        RandomLognormal
// DESCRIPTION: return a value having lognormal
//              distribution
// AUTHOR:      Tingchang YIN
// DATE:        14/04/2022
// ====================================================
__device__ __host__ float cuDFNsys::RandomLognormal(float mean,
                                                    float variance,
                                                    float min,
                                                    float max,
                                                    float rand_0_1)
{
    float mu_1 = log(mean * mean / (pow(variance + mean * mean, 0.5f)));
    float sigma_1 = pow(log(1 + ((float)variance) / (mean * mean)), 0.5f);
    //sigma_1 is input std. deviation

    float l = (log(min) - mu_1) / sigma_1;
    float u = (log(max) - mu_1) / sigma_1;

    float p_l = cuDFNsys::StandardNormalCDF(l);
    float p_u = cuDFNsys::StandardNormalCDF(u);

    float x_prime = p_l + (p_u - p_l) * rand_0_1;

    float z = cuDFNsys::StandardNormalCDFInv(x_prime);

    float value = exp(z * sigma_1 + mu_1);

    return value;
};

// ====================================================
// NAME:        StandardNormalCDF
// DESCRIPTION: return CDF of standard normal distribution (-inf, x)
// AUTHOR:      Tingchang YIN
// DATE:        14/04/2022
// ====================================================
__device__ __host__ float cuDFNsys::StandardNormalCDF(float x)
{
    float a1 = 0.398942280444;
    float a2 = 0.399903438504;
    float a3 = 5.75885480458;
    float a4 = 29.8213557808;
    float a5 = 2.62433121679;
    float a6 = 48.6959930692;
    float a7 = 5.92885724438;
    float b0 = 0.398942280385;
    float b1 = 3.8052e-08;
    float b2 = 1.00000615302;
    float b3 = 3.98064794e-04;
    float b4 = 1.98615381364;
    float b5 = 0.151679116635;
    float b6 = 5.29330324926;
    float b7 = 4.8385912808;
    float b8 = 15.1508972451;
    float b9 = 0.742380924027;
    float b10 = 30.789933034;
    float b11 = 3.99019417011;
    float cdf;
    float q;
    float y;
    //
    //  |X| <= 1.28.
    //
    if (fabs(x) <= 1.28f)
    {
        y = 0.5f * x * x;

        q = 0.5f - fabs(x) * (a1 - a2 * y / (y + a3 - a4 / (y + a5 + a6 / (y + a7))));
        //
        //  1.28 < |X| <= 12.7
        //
    }
    else if (fabs(x) <= 12.7f)
    {
        y = 0.5f * x * x;

        q = exp(-y) * b0 / (fabs(x) - b1 + b2 / (fabs(x) + b3 + b4 / (fabs(x) - b5 + b6 / (fabs(x) + b7 - b8 / (fabs(x) + b9 + b10 / (fabs(x) + b11))))));
        //
        //  12.7 < |X|
        //
    }
    else
    {
        q = 0.0f;
    }
    //
    //  Take account of negative X.
    //
    if (x < 0.0f)
    {
        cdf = q;
    }
    else
    {
        cdf = 1.0f - q;
    }

    return cdf;
};

// ====================================================
// NAME:        StandardNormalCDFInv
// DESCRIPTION: return a value that is calculated by
//              inverting the standard normal CDF.
// AUTHOR:      Tingchang YIN
// DATE:        15/04/2022
// ====================================================
__device__ __host__ float cuDFNsys::StandardNormalCDFInv(float p)
{
    float a[8] = {
        3.3871328727963666080, 1.3314166789178437745e+2,
        1.9715909503065514427e+3, 1.3731693765509461125e+4,
        4.5921953931549871457e+4, 6.7265770927008700853e+4,
        3.3430575583588128105e+4, 2.5090809287301226727e+3};
    float b[8] = {
        1.0, 4.2313330701600911252e+1,
        6.8718700749205790830e+2, 5.3941960214247511077e+3,
        2.1213794301586595867e+4, 3.9307895800092710610e+4,
        2.8729085735721942674e+4, 5.2264952788528545610e+3};
    float c[8] = {
        1.42343711074968357734, 4.63033784615654529590,
        5.76949722146069140550, 3.64784832476320460504,
        1.27045825245236838258, 2.41780725177450611770e-1,
        2.27238449892691845833e-2, 7.74545014278341407640e-4};
    float const1 = 0.180625;
    float const2 = 1.6;
    float d[8] = {
        1.0, 2.05319162663775882187,
        1.67638483018380384940, 6.89767334985100004550e-1,
        1.48103976427480074590e-1, 1.51986665636164571966e-2,
        5.47593808499534494600e-4, 1.05075007164441684324e-9};
    float e[8] = {
        6.65790464350110377720, 5.46378491116411436990,
        1.78482653991729133580, 2.96560571828504891230e-1,
        2.65321895265761230930e-2, 1.24266094738807843860e-3,
        2.71155556874348757815e-5, 2.01033439929228813265e-7};
    float f[8] = {
        1.0, 5.99832206555887937690e-1,
        1.36929880922735805310e-1, 1.48753612908506148525e-2,
        7.86869131145613259100e-4, 1.84631831751005468180e-5,
        1.42151175831644588870e-7, 2.04426310338993978564e-15};
    float q;
    float r;
    const float r8_huge = 1.0E+30;
    float split1 = 0.425;
    float split2 = 5.0;
    float value;

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

// ====================================================
// NAME:        R8PolyValueHorner
// DESCRIPTION: A stable algorithm to calculate
//              polynomial with 8th degrees
// AUTHOR:      Tingchang YIN
// DATE:        15/04/2022
// ====================================================
__device__ __host__ float cuDFNsys::R8PolyValueHorner(int m, float c[], float x)
{
    int i;
    float value;

    value = c[m];

    for (i = m - 1; 0 <= i; i--)
    {
        value = value * x + c[i];
    }

    return value;
};