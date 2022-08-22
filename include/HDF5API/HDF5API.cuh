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
    vector<double> ReadDataset(const string &filename,
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
                             const string &datasetname)
    {
        try
        {
            H5::Exception::dontPrint();
            H5File file1(filename, H5F_ACC_RDONLY);
            file1.close();
        }
        catch (...)
        {
            string AS = "File '" + filename + "' does not exist!\n";
            throw ExceptionsPause(AS);
        };

        H5File file(filename, H5F_ACC_RDONLY);
        DataSet dataset;
        Group group;
        if (groupname != "N")
        {
            Group group = file.openGroup(groupname);
            dataset = group.openDataSet(datasetname);
        }
        else
            dataset = file.openDataSet(datasetname);

        DataSpace filespace = dataset.getSpace();
        H5::StrType datatype = dataset.getStrType();

        std::string data;

        dataset.read(data, datatype, filespace);

        if (groupname != "N")
            group.close();
        file.close();

        return data;
    };
};
}; // namespace cuDFNsys