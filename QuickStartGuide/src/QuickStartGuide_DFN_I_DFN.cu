#include "cuDFNsys.cuh"
int main()
{
    time_t t;
    time(&t);
    cuDFNsys::DFN<double> my_dfn;
    my_dfn.NumFractures = {70, 80};
    my_dfn.Kappa = {20, 10};
    my_dfn.MeanOrientationOfFisherDistribution = {make_double3(0., 0., 1.),
                                                  make_double3(1., 0., 0.)};
    my_dfn.DomainSizeX = 30;
    my_dfn.DomainDimensionRatio = make_double3(1., 1., 2.);
    my_dfn.Beta = {0.2, 0.3};
    my_dfn.Gamma = {2.0e-5, 3.0e-6};
    my_dfn.ModeOfSizeDistribution = {0, 1};
    my_dfn.SizeDistributionParameters = {make_double4(1.5, 1., 15., 0.),
                                         make_double4(8.5, 5.5, 1., 15.)};
    my_dfn.PercoDir = 2;
    my_dfn.RandomSeed = (unsigned long)t;
    my_dfn.FractureGeneration();
    my_dfn.IdentifyIntersectionsClusters(true);
    my_dfn.Visualization("DFN_VISUAL_I", "DFN_VISUAL_I", "DFN_VISUAL_I", true,
                         true, true, true);
    my_dfn.StoreInH5("Class_DFN");
    return 0;
};