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
// NAME:        DispersionAtOneDensityValue.cu
// DESCRIPTION: Dispersion in a DFN with a specific percolation parameter value
// AUTHOR:      Tingchang YIN
// DATE:        24/03/2023
// ====================================================

#include "cuDFNsys.cuh"
#include <fstream>
#include <iostream>
#include <unistd.h>
#ifdef USE_DOUBLES
typedef double _DataType_;
#else
typedef float _DataType_;
#endif

int main(int argc, char *argv[])
{
    for (uint i = 0; i < 2; ++i)
    {
        string path2 = "a" + std::to_string(i);
        string command1 = "mkdir -p " + path2;
        system(command1.c_str());

        std::ofstream fs("./SimulationFailed" + std::to_string(i) + ".txt");
        fs.close();

        
    }

    return 0;
};
