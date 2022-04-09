#include "HDF5API/HDF5API.cuh"

// ====================================================
// NAME:        NewFile
// DESCRIPTION: create a new h5 file
// AUTHOR:      Tingchang YIN
// DATE:        09/04/2022
// ====================================================
void cuDFNsys::HDF5API::NewFile(const string &name)
{
    try
    {
        H5::Exception::dontPrint();
        H5File file(name, H5F_ACC_TRUNC);
        file.close();
    }
    catch (...)
    {
        cout << "\033[31mA file with same name does exist and is occupied by anthor program now!\033[0m\n";
        cout << "\033[31mSo, I will delete this file and create a new one!\033[0m\n";
        string name1 = name;
        std::remove(name1.c_str());
        H5File file(name, H5F_ACC_TRUNC);
        file.close();
    }
}; // NewFile

// ====================================================
// NAME:        AddDataset
// DESCRIPTION: add dataset to a h5 file
// AUTHOR:      Tingchang YIN
// DATE:        09/04/2022
// ====================================================
template <class T>
void cuDFNsys::HDF5API::AddDataset(const string &filename,
                                   const string &groupname,
                                   const string &datasetname,
                                   const T *data, // the T should be column-major
                                   const uint2 &dim)
{
    H5File file(filename, H5F_ACC_RDWR);

    hsize_t *dims; // dataset dimensions for each rank
    int rank_ = 1;
    if (dim.y > 1)
        rank_ = 2;

    dims = new hsize_t[rank_];
    dims[0] = dim.x;
    if (rank_ == 2)
        dims[1] = dim.y;
    DataSpace dataspace(rank_, dims);
    delete[] dims;
    dims = NULL;

    if (groupname != "N")
    {
        Group group;

        H5::Exception::dontPrint();
        try
        {
            //cout << "try to open a group!\n";
            group = file.openGroup(groupname);
            //cout << "opened group!\n";
        }
        catch (...)
        {
            //cout << "no this group! create a new group!\n";
            group = file.createGroup(groupname);
            //cout << "created group!\n";
        }

        DataSet dataset =
            group.createDataSet(datasetname, PredType::NATIVE_DOUBLE, dataspace);

        double *buffer = new double[dim.x * dim.y]();
        if (buffer == NULL)
        {
            string AS = "Alloc error in HDF5API::AddDataset!\n";
            group.close();
            file.close();
            throw ExceptionsPause(AS);
        }

        for (size_t i = 0; i < dim.x * dim.y; ++i)
            buffer[i] = data[i];

        dataset.write(buffer, PredType::NATIVE_DOUBLE);

        delete[] buffer;
        buffer = NULL;

        group.close();
    }
    else
    {
        DataSet dataset =
            file.createDataSet(datasetname, PredType::NATIVE_DOUBLE, dataspace);

        double *buffer = new double[dim.x * dim.y]();
        if (buffer == NULL)
        {
            string AS = "Alloc error in HDF5API::AddDataset!\n";
            file.close();
            throw ExceptionsPause(AS);
        }

        for (size_t i = 0; i < dim.x * dim.y; ++i)
            buffer[i] = data[i];

        dataset.write(buffer, PredType::NATIVE_DOUBLE);

        delete[] buffer;
        buffer = NULL;
    }

    file.close();
}; // AddDataset
template void cuDFNsys::HDF5API::AddDataset(const string &filename, const string &groupname, const string &datasetname, const int *data, const uint2 &dim);
template void cuDFNsys::HDF5API::AddDataset(const string &filename, const string &groupname, const string &datasetname, const double *data, const uint2 &dim);
template void cuDFNsys::HDF5API::AddDataset(const string &filename, const string &groupname, const string &datasetname, const float *data, const uint2 &dim);
template void cuDFNsys::HDF5API::AddDataset(const string &filename, const string &groupname, const string &datasetname, const size_t *data, const uint2 &dim);
template void cuDFNsys::HDF5API::AddDataset(const string &filename, const string &groupname, const string &datasetname, const uint *data, const uint2 &dim);

// ====================================================
// NAME:        AddDataset
// DESCRIPTION: add a number of datasets with a group
//              name to a h5 file
// AUTHOR:      Tingchang YIN
// DATE:        09/04/2022
// ====================================================
template <class T>
void cuDFNsys::HDF5API::AddDatasetsWithOneGroup(const string &filename,
                                                const string &groupname,
                                                const vector<string> &datasetname,
                                                const vector<T *> data,
                                                const vector<uint2> &dim)
{
    H5File file(filename, H5F_ACC_RDWR);

    Group group;

    H5::Exception::dontPrint();
    try
    {
        //cout << "try to open a group!\n";
        group = file.openGroup(groupname);
        //cout << "opened group!\n";
    }
    catch (...)
    {
        //cout << "no this group! create a new group!\n";
        group = file.createGroup(groupname);
        //cout << "created group!\n";
    };

    for (int i = 0; i < datasetname.size(); ++i)
    {
        hsize_t *dims; // dataset dimensions for each rank
        int rank_ = 1;
        if (dim[i].y > 1)
            rank_ = 2;

        dims = new hsize_t[rank_];
        dims[0] = dim[i].x;
        if (rank_ == 2)
            dims[1] = dim[i].y;
        DataSpace dataspace(rank_, dims);
        delete[] dims;
        dims = NULL;

        DataSet dataset =
            group.createDataSet(datasetname[i], PredType::NATIVE_DOUBLE, dataspace);

        double *buffer = new double[dim[i].x * dim[i].y]();
        if (buffer == NULL)
        {
            string AS = "Alloc error in HDF5API::AddDatasetsWithOneGroup!\n";
            group.close();
            file.close();
            throw ExceptionsPause(AS);
        }

        for (size_t j = 0; j < dim[i].x * dim[i].y; ++j)
            buffer[j] = data[i][j];

        dataset.write(buffer, PredType::NATIVE_DOUBLE);

        delete[] buffer;
        buffer = NULL;
    }

    group.close();
    file.close();
}; // AddDatasetsWithOneGroup
template void cuDFNsys::HDF5API::AddDatasetsWithOneGroup(const string &filename, const string &groupname, const vector<string> &datasetname, const vector<int *> data, const vector<uint2> &dim);
template void cuDFNsys::HDF5API::AddDatasetsWithOneGroup(const string &filename, const string &groupname, const vector<string> &datasetname, const vector<double *> data, const vector<uint2> &dim);
template void cuDFNsys::HDF5API::AddDatasetsWithOneGroup(const string &filename, const string &groupname, const vector<string> &datasetname, const vector<float *> data, const vector<uint2> &dim);
template void cuDFNsys::HDF5API::AddDatasetsWithOneGroup(const string &filename, const string &groupname, const vector<string> &datasetname, const vector<size_t *> data, const vector<uint2> &dim);
template void cuDFNsys::HDF5API::AddDatasetsWithOneGroup(const string &filename, const string &groupname, const vector<string> &datasetname, const vector<uint *> data, const vector<uint2> &dim);

// ====================================================
// NAME:        OverWrite
// DESCRIPTION: OverWrite dataset to a h5
// AUTHOR:      Tingchang YIN
// DATE:        09/04/2022
// ====================================================
template <class T>
void cuDFNsys::HDF5API::OverWrite(const string &filename,
                                  const string &groupname,
                                  const string &datasetname,
                                  const T *data,
                                  const uint2 &dim)
{
    H5File file(filename, H5F_ACC_RDWR); //The hdf5 c++ object.

    string channelName;
    if (groupname != "N")
        channelName = "/" + groupname + "/" + datasetname;
    else
        channelName = "/" + datasetname;

    int result = H5Ldelete(file.getId(), channelName.data(), H5P_DEFAULT);
    result++;
    file.close();

    this->AddDataset(filename, groupname,
                     datasetname, data, dim);
}; // OverWrite
template void cuDFNsys::HDF5API::OverWrite(const string &filename, const string &groupname, const string &datasetname, const int *data, const uint2 &dim);
template void cuDFNsys::HDF5API::OverWrite(const string &filename, const string &groupname, const string &datasetname, const double *data, const uint2 &dim);
template void cuDFNsys::HDF5API::OverWrite(const string &filename, const string &groupname, const string &datasetname, const float *data, const uint2 &dim);
template void cuDFNsys::HDF5API::OverWrite(const string &filename, const string &groupname, const string &datasetname, const size_t *data, const uint2 &dim);
template void cuDFNsys::HDF5API::OverWrite(const string &filename, const string &groupname, const string &datasetname, const uint *data, const uint2 &dim);

// ====================================================
// NAME:        ReadDataset
// DESCRIPTION: Read dataset in a h5
// AUTHOR:      Tingchang YIN
// DATE:        09/04/2022
// ====================================================
vector<double> cuDFNsys::HDF5API::ReadDataset(const string &filename,
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
    int rank = filespace.getSimpleExtentNdims();

    hsize_t dims[rank];

    rank = filespace.getSimpleExtentDims(dims);

    DataSpace myspace(rank, dims);

    int NUM_size = 1;
    for (int i = 0; i < rank; ++i)
        NUM_size *= dims[i];

    double *buffer = new double[NUM_size]();
    if (buffer == NULL)
    {
        string AS = "Alloc error in HDF5API::ReadDataset!\n";
        throw ExceptionsPause(AS);
    }

    dataset.read(buffer, PredType::NATIVE_DOUBLE, myspace, filespace);

    vector<double> AK(buffer, buffer + NUM_size);

    delete[] buffer;
    buffer = NULL;
    if (groupname != "N")
        group.close();
    file.close();

    return AK;
}; // ReadDataset

// ====================================================
// NAME:        IfH5FileExist
// DESCRIPTION: check if a h5 exists?
// AUTHOR:      Tingchang YIN
// DATE:        09/04/2022
// ====================================================
bool cuDFNsys::HDF5API::IfH5FileExist(const string &filename)
{
    try
    {
        H5::Exception::dontPrint();
        H5File file(filename, H5F_ACC_RDONLY);
        file.close();
        return true;
    }
    catch (...)
    {
        return false;
    };
    return false;
}; //IfH5FileExist