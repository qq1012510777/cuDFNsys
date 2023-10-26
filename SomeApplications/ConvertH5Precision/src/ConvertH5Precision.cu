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
// NAME:        A test case
// DESCRIPTION: Call cuDFNsys functions to do simulation and test.
// AUTHOR:      Tingchang YIN
// DATE:        30/06/2022
// ====================================================

#include "cuDFNsys.cuh"
#include <unistd.h>

#ifdef USE_DOUBLES
typedef double _DataType_;
#else
typedef float _DataType_;
#endif

int main(int argc, char *argv[])
{

    try
    {
        cuDFNsys::HDF5API h5g;
        string filename = argv[1];
        string datasetname = argv[2];
        vector<float> dataFDS = h5g.ReadDataset<float>(filename,
                                                       "N", datasetname);

        string newfile = "newH5.h5";
        h5g.NewFile(newfile);

        //--------------------
        H5File file(newfile, H5F_ACC_RDWR);

        hsize_t *dims; // dataset dimensions for each rank
        dims = new hsize_t[3];

        dims[0] = atoi(argv[3]);
        dims[1] = atoi(argv[4]);
        dims[2] = atoi(argv[5]);

        DataSpace dataspace(3, dims);
        delete[] dims;
        dims = NULL;

        auto datatype = PredType::NATIVE_FLOAT;

        DataSet dataset =
            file.createDataSet(datasetname, datatype, dataspace);

        dataset.write(dataFDS.data(), datatype);

        file.close();
    }
    catch (cuDFNsys::ExceptionsIgnore &e)
    {
        cout << e.what() << endl;
    }
    catch (cuDFNsys::ExceptionsPause &e)
    {
        cout << e.what() << endl;
    }
    catch (...)
    {
        throw;
    };
    return 0;
};
