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