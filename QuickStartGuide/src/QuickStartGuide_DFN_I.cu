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
// NAME:        A quickstart example to generate DFNs
// DESCRIPTION: Call cuDFNsys functions to do simulation.
// AUTHOR:      Tingchang YIN
// DATE:        13/10/2023
// ====================================================

#include "cuDFNsys.cuh"
#include <fstream>
#include <iostream>
#include <limits.h>
#include <unistd.h>

int main(int argc, char *argv[])
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

    try
    {
        cuDFNsys::DFN<double> my_dfn222222;
        my_dfn222222.LoadClassFromH5("Class_DFN");
        my_dfn222222.Visualization("DFN_VISUAL_II", "DFN_VISUAL_II",
                                   "DFN_VISUAL_II", true, true, true, true);
    }
    catch (cuDFNsys::ExceptionsPause &e)
    {
        cout << e.what() << endl;
        throw;
    }
    catch (H5::Exception &e)
    {
        cout << "H5::Exception\n";
        std::cout << e.getDetailMsg() << std::endl;
        throw;
    }

    return 0;
};