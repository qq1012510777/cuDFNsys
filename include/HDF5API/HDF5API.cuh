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

///////////////////////////////////////////////////////////////////
// NAME:              HDF5API.h
//
// PURPOSE:           an h5df api to process data file
//
// FUNCTIONS/OBJECTS: HDF5API
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../Exceptions/Exceptions.cuh"
#include "H5Cpp.h"
#include "hdf5.h"
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using namespace H5;

namespace cuDFNsys
{
class HDF5API
{
public:
    //constructor
    HDF5API(){};
    // new file
    void NewFile(const string &name);
    // add a new dataset with a new group, if group = "N", no group
    template <class T>
    void AddDataset(const string &filename,
                    const string &groupname,
                    const string &datasetname,
                    const T *data, // the T should be column-major
                    const uint2 &dim);
    // add a number of datasets to a h5, datasets are within a group
    template <class T>
    void AddDatasetsWithOneGroup(const string &filename,
                                 const string &groupname,
                                 const vector<string> &datasetname,
                                 const vector<T *> data,
                                 const vector<uint2> &dim);
    // overwrite a dataset in a h5
    template <class T>
    void OverWrite(const string &filename,
                   const string &groupname,
                   const string &datasetname,
                   const T *data,
                   const uint2 &dim);
    // read dataset in a h5; the function returns a vector that is column major
    template <class T>
    vector<T> ReadDataset(const string &filename,
                          const string &groupname,
                          const string &datasetname);
    //check if a h5 exits?
    bool IfH5FileExist(const string &filename);
    // add string data
    void AddDatasetString(const string &filename,
                          const string &groupname,
                          const string &datasetname,
                          const string &sdata);
    // read string
    string ReadDatasetString(const string &filename,
                             const string &groupname,
                             const string &datasetname);
};
}; // namespace cuDFNsys