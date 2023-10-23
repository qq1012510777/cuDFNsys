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

#include "InputObjectData/InputObjectData.cuh"

// ====================================================
// NAME:        InputFractures
// DESCRIPTION: InputFractures
// AUTHOR:      Tingchang YIN
// DATE:        02/08/2022
// ====================================================
template <typename T>
void cuDFNsys::InputObjectData<T>::InputFractures(const string &filename,
                                                  thrust::host_vector<cuDFNsys::Fracture<T>> &Frac_verts_host,
                                                  T &L, double3 &DomainDimensionRatio)
{
    //

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

    uint Dsize = 0;

    H5File file(filename, H5F_ACC_RDONLY);
    DataSet dataset;
    dataset = file.openDataSet("NumFractures");

    DataSpace filespace = dataset.getSpace();
    int rank = filespace.getSimpleExtentNdims();

    hsize_t dims[rank];

    rank = filespace.getSimpleExtentDims(dims);

    DataSpace myspace(rank, dims);

    int NUM_size = 1;
    for (int i = 0; i < rank; ++i)
        NUM_size *= dims[i];

    if (NUM_size != 1)
    {
        string AS = "Number of fractures should be equal to one!\n";
        throw ExceptionsPause(AS);
    }

    double *buffer = new double[NUM_size]();
    if (buffer == NULL)
    {
        string AS = "Alloc error in InputFractures::InputFractures!\n";
        throw ExceptionsPause(AS);
    }

    dataset.read(buffer, PredType::NATIVE_DOUBLE, myspace, filespace);

    Dsize = (uint)buffer[0];

    delete[] buffer;
    buffer = NULL;

    //--------------------------- finish getting the number of fractures
    //--------------------------- finish getting the number of fractures
    //--------------------------- finish getting the number of fractures

    DataSet dataset2 = file.openDataSet("L");
    DataSpace filespace2 = dataset2.getSpace();
    int rank2 = filespace2.getSimpleExtentNdims();

    hsize_t dims2[rank2];

    rank2 = filespace2.getSimpleExtentDims(dims2);

    DataSpace myspace2(rank2, dims2);

    int NUM_size2 = 1;
    for (int i = 0; i < rank2; ++i)
        NUM_size2 *= dims2[i];

    if (NUM_size2 != 1)
    {
        string AS = "dimension should be equal to one!\n";
        throw ExceptionsPause(AS);
    }

    double *buffer2 = new double[NUM_size2]();
    if (buffer2 == NULL)
    {
        string AS = "Alloc error in InputFractures::InputFractures!\n";
        throw ExceptionsPause(AS);
    }

    dataset2.read(buffer2, PredType::NATIVE_DOUBLE, myspace2, filespace2);

    L = (uint)buffer2[0];

    delete[] buffer2;
    buffer2 = NULL;

    //--------------------------- finish getting model size
    //--------------------------- finish getting model size
    //--------------------------- finish getting model size
    try
    {
        cuDFNsys::HDF5API h5g;

        vector<double> Di_ = h5g.ReadDataset<double>(filename,
                                                     "N", "DomainDimensionRatio");
        DomainDimensionRatio.x = Di_[0];
        DomainDimensionRatio.y = Di_[1];
        DomainDimensionRatio.z = Di_[2];
    }
    catch (...)
    {
        DomainDimensionRatio = make_double3(1, 1, 1);
    }
    //----------------------------------

    Frac_verts_host.resize(Dsize);

    for (uint i = 0; i < Dsize; ++i)
    {
        string groupname = "Fracture_" + to_string(i + 1);
        Group group = file.openGroup(groupname);
        //cout << "Dsize = " << Dsize << ", " << groupname << endl;
        vector<string> dataset_string_vec = {"Conductivity",
                                             "Verts3D",
                                             "Center",
                                             "Verts3DTruncated",
                                             "NumVertsTruncated",
                                             "Radius",
                                             "ConnectModelSurf",
                                             "NormalVec"};
        T Conductivity_[1];
        T Verts3D_[12];
        T Center_[3];
        T Verts3DTruncated_[24];
        T NumVertsTruncated_[1];
        T Radius_[1];
        T ConnectModelSurf_[6];
        T NormalVec_[3];

        vector<T *> data = {Conductivity_,
                            Verts3D_,
                            Center_,
                            Verts3DTruncated_,
                            NumVertsTruncated_,
                            Radius_,
                            ConnectModelSurf_,
                            NormalVec_};

        for (uint j = 0; j < dataset_string_vec.size(); ++j)
        {
            DataSet dataset_1 = group.openDataSet(dataset_string_vec[j]);

            DataSpace filespace_1 = dataset_1.getSpace();
            int rank_1 = filespace_1.getSimpleExtentNdims();

            hsize_t dims_1[rank_1];

            rank_1 = filespace_1.getSimpleExtentDims(dims_1);

            DataSpace myspace_1(rank_1, dims_1);

            int NUM_size = 1;
            for (int ik = 0; ik < rank_1; ++ik)
                NUM_size *= dims_1[ik];

            double *buffer = new double[NUM_size]();
            if (buffer == NULL)
            {
                string AS = "Alloc error in InputFractures::InputFractures!\n";
                throw ExceptionsPause(AS);
            }

            dataset_1.read(buffer, PredType::NATIVE_DOUBLE, myspace_1, filespace_1);

            for (uint ik = 0; ik < NUM_size; ++ik)
                data[j][ik] = buffer[ik];

            delete[] buffer;
            buffer = NULL;
        }
        group.close();

        Frac_verts_host[i].Conductivity = Conductivity_[0];

        Frac_verts_host[i].Verts3D[0] = cuDFNsys::MakeVector3(Verts3D_[0], Verts3D_[4], Verts3D_[8]);
        Frac_verts_host[i].Verts3D[1] = cuDFNsys::MakeVector3(Verts3D_[1], Verts3D_[5], Verts3D_[9]);
        Frac_verts_host[i].Verts3D[2] = cuDFNsys::MakeVector3(Verts3D_[2], Verts3D_[6], Verts3D_[10]);
        Frac_verts_host[i].Verts3D[3] = cuDFNsys::MakeVector3(Verts3D_[3], Verts3D_[7], Verts3D_[11]);

        Frac_verts_host[i].Center = cuDFNsys::MakeVector3(Center_[0], Center_[1], Center_[2]);

        Frac_verts_host[i].Verts3DTruncated[0] = cuDFNsys::MakeVector3(Verts3DTruncated_[0], Verts3DTruncated_[8], Verts3DTruncated_[16]);
        Frac_verts_host[i].Verts3DTruncated[1] = cuDFNsys::MakeVector3(Verts3DTruncated_[1], Verts3DTruncated_[9], Verts3DTruncated_[17]);
        Frac_verts_host[i].Verts3DTruncated[2] = cuDFNsys::MakeVector3(Verts3DTruncated_[2], Verts3DTruncated_[10], Verts3DTruncated_[18]);
        Frac_verts_host[i].Verts3DTruncated[3] = cuDFNsys::MakeVector3(Verts3DTruncated_[3], Verts3DTruncated_[11], Verts3DTruncated_[19]);
        Frac_verts_host[i].Verts3DTruncated[4] = cuDFNsys::MakeVector3(Verts3DTruncated_[4], Verts3DTruncated_[12], Verts3DTruncated_[20]);
        Frac_verts_host[i].Verts3DTruncated[5] = cuDFNsys::MakeVector3(Verts3DTruncated_[5], Verts3DTruncated_[13], Verts3DTruncated_[21]);
        Frac_verts_host[i].Verts3DTruncated[6] = cuDFNsys::MakeVector3(Verts3DTruncated_[6], Verts3DTruncated_[14], Verts3DTruncated_[22]);
        Frac_verts_host[i].Verts3DTruncated[7] = cuDFNsys::MakeVector3(Verts3DTruncated_[7], Verts3DTruncated_[15], Verts3DTruncated_[23]);

        Frac_verts_host[i].NumVertsTruncated = NumVertsTruncated_[0];

        Frac_verts_host[i].Radius = Radius_[0];

        Frac_verts_host[i].ConnectModelSurf[0] = ConnectModelSurf_[0];
        Frac_verts_host[i].ConnectModelSurf[1] = ConnectModelSurf_[1];
        Frac_verts_host[i].ConnectModelSurf[2] = ConnectModelSurf_[2];
        Frac_verts_host[i].ConnectModelSurf[3] = ConnectModelSurf_[3];
        Frac_verts_host[i].ConnectModelSurf[4] = ConnectModelSurf_[4];
        Frac_verts_host[i].ConnectModelSurf[5] = ConnectModelSurf_[5];

        Frac_verts_host[i].NormalVec = cuDFNsys::MakeVector3(NormalVec_[0], NormalVec_[1], NormalVec_[2]);
    }

    file.close();

}; // InputFractures
template void cuDFNsys::InputObjectData<double>::InputFractures(const string &filename,
                                                                thrust::host_vector<cuDFNsys::Fracture<double>> &Frac_verts_host,
                                                                double &L, double3 &DomainDimensionRatio);
template void cuDFNsys::InputObjectData<float>::InputFractures(const string &filename,
                                                               thrust::host_vector<cuDFNsys::Fracture<float>> &Frac_verts_host,
                                                               float &L, double3 &DomainDimensionRatio);

// ====================================================
// NAME:        InputMesh
// DESCRIPTION: InputMesh
// AUTHOR:      Tingchang YIN
// DATE:        02/08/2022
// ====================================================
template <typename T>
void cuDFNsys::InputObjectData<T>::InputMesh(const string &filename,
                                             cuDFNsys::Mesh<T> &mesh,
                                             std::vector<size_t> *Fracs_percol)
{
    //
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

    HDF5API h5g;
    std::vector<double> TMP_U = h5g.ReadDataset<double>(filename, "N", "MeanGridSize");
    mesh.MeanGridSize = TMP_U[0];

    H5File file(filename, H5F_ACC_RDONLY);

    vector<string> datasetname = {
        "NumFractures",
        "NumNodes",
        "NumElements",
        "NumInletEdges",
        "NumOutletEdges"};

    uint NumFractures = 0,
         NumNodes = 0,
         NumElements = 0,
         NumInletEdges = 0,
         NumOutletEdges = 0;

    vector<uint *> data = {
        &NumFractures,
        &NumNodes,
        &NumElements,
        &NumInletEdges,
        &NumOutletEdges};

    for (uint i = 0; i < datasetname.size(); ++i)
    {

        DataSet dataset;
        dataset = file.openDataSet(datasetname[i]);

        DataSpace filespace = dataset.getSpace();
        int rank = filespace.getSimpleExtentNdims();

        hsize_t dims[rank];

        rank = filespace.getSimpleExtentDims(dims);

        DataSpace myspace(rank, dims);

        int NUM_size = 1;
        for (int j = 0; j < rank; ++j)
            NUM_size *= dims[j];

        if (NUM_size != 1)
        {
            string AS = "Number of fractures should be equal to one!\n";
            throw ExceptionsPause(AS);
        }

        double *buffer = new double[NUM_size]();
        if (buffer == NULL)
        {
            string AS = "Alloc error in InputFractures::InputFractures!\n";
            throw ExceptionsPause(AS);
        }

        dataset.read(buffer, PredType::NATIVE_DOUBLE, myspace, filespace);

        *(data[i]) = (uint)buffer[0];
        //cout << *(data[i]) << endl;
        delete[] buffer;
        buffer = NULL;
    }

    /// ------------------- group_mesh----------------
    /// ------------------- group_mesh----------------
    /// ------------------- group_mesh----------------
    string groupname = "group_mesh";
    Group group = file.openGroup(groupname);

    vector<string> datasetname_x = {
        "MeshSuccess",
        "FracID",
        "Coordinate3D",
        "Element2D",
        "NumElementFrac",
        "Element3D",
        "Coordinate2D",
        "ElementFracTag",
        "EdgeAttri",
        "InletEdgeNOLen",
        "OutletEdgeNOLen",
        "NumInteriorEdges",
        "NumNeumannEdges",
        "Dir"};

    T MeshSuccess[1];
    T *FracID = new T[NumFractures];
    T *Coordinate3D = new T[3 * NumNodes];
    T *Element2D = new T[3 * NumElements];
    T *NumElementFrac = new T[NumFractures];
    T *Element3D = new T[3 * NumElements];
    T *Coordinate2D = new T[NumElements * 6];
    T *ElementFracTag = new T[NumElements];
    T *EdgeAttri = new T[6 * NumElements];
    T *InletEdgeNOLen = new T[NumInletEdges * 2];
    T *OutletEdgeNOLen = new T[NumOutletEdges * 2];
    T NumInteriorEdges[1];
    T NumNeumannEdges[1];
    T Dir[1];

    vector<T *> data_x = {
        MeshSuccess,
        FracID,
        Coordinate3D,
        Element2D,
        NumElementFrac,
        Element3D,
        Coordinate2D,
        ElementFracTag,
        EdgeAttri,
        InletEdgeNOLen,
        OutletEdgeNOLen,
        NumInteriorEdges,
        NumNeumannEdges,
        Dir};

    for (uint i = 0; i < datasetname_x.size(); ++i)
    {
        DataSet dataset_1 = group.openDataSet(datasetname_x[i]);

        DataSpace filespace_1 = dataset_1.getSpace();
        int rank_1 = filespace_1.getSimpleExtentNdims();

        hsize_t dims_1[rank_1];

        rank_1 = filespace_1.getSimpleExtentDims(dims_1);

        DataSpace myspace_1(rank_1, dims_1);

        int NUM_size = 1;
        for (int ik = 0; ik < rank_1; ++ik)
            NUM_size *= dims_1[ik];

        double *buffer = new double[NUM_size]();
        if (buffer == NULL)
        {
            string AS = "Alloc error in InputFractures::InputFractures!\n";
            throw ExceptionsPause(AS);
        }

        dataset_1.read(buffer, PredType::NATIVE_DOUBLE, myspace_1, filespace_1);

        for (uint ik = 0; ik < NUM_size; ++ik)
            data_x[i][ik] = buffer[ik];

        delete[] buffer;
        buffer = NULL;
    }

    if (MeshSuccess[0] == 1)
        mesh.MeshSuccess = true;
    else
        mesh.MeshSuccess = false;

    mesh.FracID = Fracs_percol;

    mesh.Coordinate3D.resize(NumNodes);
    for (uint i = 0; i < 3; ++i)
        for (uint j = 0; j < NumNodes; ++j)
        {
            T *KL = &(mesh.Coordinate3D[j].x);
            KL[i] = Coordinate3D[j + i * NumNodes];
        }

    mesh.Element2D.resize(NumFractures);
    for (uint i = 0; i < NumFractures; ++i)
        mesh.Element2D[i].resize(NumElementFrac[i]);

    uint iio = 0;
    for (uint i = 0; i < 3; ++i)
        for (uint j = 0; j < NumFractures; ++j)
            for (uint k = 0; k < mesh.Element2D[j].size(); ++k)
            {
                uint *KL = &(mesh.Element2D[j][k].x);
                KL[i] = (uint)Element2D[iio];
                iio++;
            };

    mesh.Element3D.resize(NumElements);
    for (uint i = 0; i < 3; ++i)
        for (uint j = 0; j < NumElements; ++j)
        {
            uint *KL = &(mesh.Element3D[j].x);
            KL[i] = (uint)Element3D[i * NumElements + j];
        }

    mesh.Coordinate2D.resize(NumElements);
    for (uint i = 0; i < 6; ++i)
        for (uint j = 0; j < NumElements; ++j)
        {
            T *KL;

            if (i < 3)
                KL = &(mesh.Coordinate2D[j].x[0]);
            else
                KL = &(mesh.Coordinate2D[j].y[0]);

            KL[i % 3] = Coordinate2D[i * NumElements + j];
        }

    mesh.ElementFracTag.resize(NumElements);
    for (uint i = 0; i < NumElements; ++i)
        mesh.ElementFracTag[i] = (uint)ElementFracTag[i];

    mesh.EdgeAttri.resize(NumElements);
    for (uint i = 0; i < 6; ++i)
        for (uint j = 0; j < NumElements; ++j)
        {
            int *KL;

            if (i < 3)
                KL = &(mesh.EdgeAttri[j].e[0]);
            else
                KL = &(mesh.EdgeAttri[j].no[0]);

            KL[i % 3] = (int)EdgeAttri[i * NumElements + j];
        }

    mesh.InletEdgeNOLen.resize(NumInletEdges);
    for (uint i = 0; i < 2; ++i)
        for (uint j = 0; j < NumInletEdges; ++j)
        {
            T *KL = &(mesh.InletEdgeNOLen[j].x);
            KL[i] = (T)InletEdgeNOLen[i * NumInletEdges + j];
        }

    mesh.OutletEdgeNOLen.resize(NumOutletEdges);
    for (uint i = 0; i < 2; ++i)
        for (uint j = 0; j < NumOutletEdges; ++j)
        {
            T *KL = &(mesh.OutletEdgeNOLen[j].x);
            KL[i] = (T)OutletEdgeNOLen[i * NumOutletEdges + j];
        }

    mesh.NumInteriorEdges = NumInteriorEdges[0];
    mesh.NumNeumannEdges = NumNeumannEdges[0];

    mesh.NumInletEdges = NumInletEdges;
    mesh.NumOutletEdges = NumOutletEdges;

    group.close();
    file.close();

    delete[] OutletEdgeNOLen;
    OutletEdgeNOLen = NULL;

    delete[] InletEdgeNOLen;
    InletEdgeNOLen = NULL;

    delete[] EdgeAttri;
    EdgeAttri = NULL;

    delete[] ElementFracTag;
    ElementFracTag = NULL;

    delete[] Coordinate2D;
    Coordinate2D = NULL;

    delete[] Element3D;
    Element3D = NULL;

    delete[] NumElementFrac;
    NumElementFrac = NULL;

    delete[] Element2D;
    Element2D = NULL;

    delete[] Coordinate3D;
    Coordinate3D = NULL;

    delete[] FracID;
    FracID = NULL;
}; // InputMesh
template void cuDFNsys::InputObjectData<double>::InputMesh(const string &filename,
                                                           cuDFNsys::Mesh<double> &mesh,
                                                           std::vector<size_t> *Fracs_percol);
template void cuDFNsys::InputObjectData<float>::InputMesh(const string &filename,
                                                          cuDFNsys::Mesh<float> &mesh,
                                                          std::vector<size_t> *Fracs_percol);

// ====================================================
// NAME:        InputMHFEM
// DESCRIPTION: InputMHFEM
// AUTHOR:      Tingchang YIN
// DATE:        29/12/2022
// ====================================================
template <typename T>
void cuDFNsys::InputObjectData<T>::InputMHFEM(const string &filename,
                                              cuDFNsys::MHFEM<T> &MHFEM)
{
    //
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

    HDF5API h5g;

    vector<T> value_ = h5g.ReadDataset<T>(filename, "N", "QIn");
    MHFEM.QIn = value_[0];

    value_ = h5g.ReadDataset<T>(filename, "N", "QIn");
    MHFEM.QIn = value_[0];

    value_ = h5g.ReadDataset<T>(filename, "N", "QOut");
    MHFEM.QOut = value_[0];

    value_ = h5g.ReadDataset<T>(filename, "N", "InletLength");
    MHFEM.InletLength = value_[0];

    value_ = h5g.ReadDataset<T>(filename, "N", "OutletLength");
    MHFEM.OutletLength = value_[0];

    value_ = h5g.ReadDataset<T>(filename, "N", "QError");
    MHFEM.QError = value_[0];

    value_ = h5g.ReadDataset<T>(filename, "N", "Permeability");
    MHFEM.Permeability = value_[0];

    value_ = h5g.ReadDataset<T>(filename, "N", "InletP");
    MHFEM.InletP = value_[0];

    value_ = h5g.ReadDataset<T>(filename, "N", "OutletP");
    MHFEM.OutletP = value_[0];

    vector<double> value_2 = h5g.ReadDataset<double>(filename, "N", "TripletTime");
    MHFEM.TripletTime = value_2[0];

    value_2 = h5g.ReadDataset<double>(filename, "N", "MeanVelocity");
    MHFEM.MeanVelocity = value_2[0];

    value_2 = h5g.ReadDataset<double>(filename, "N", "MaxVelocity");
    MHFEM.MaxVelocity = value_2[0];

    vector<int> value_3 = h5g.ReadDataset<int>(filename, "N", "Dir");
    MHFEM.Dir = value_3[0];

    value_2 = h5g.ReadDataset<double>(filename, "N", "PressureInteriorEdge");
    MHFEM.PressureInteriorEdge = Eigen::Map<Eigen::MatrixXd>(value_2.data(),
                                                             value_2.size(), 1);

    value_2 = h5g.ReadDataset<double>(filename, "N", "PressureEles");
    MHFEM.PressureEles = Eigen::Map<Eigen::MatrixXd>(value_2.data(),
                                                     value_2.size(), 1);

    value_2 = h5g.ReadDataset<double>(filename, "N", "VelocityNormalScalarSepEdges");
    MHFEM.VelocityNormalScalarSepEdges = Eigen::Map<Eigen::MatrixXd>(value_2.data(),
                                                                     value_2.size(), 1);
}; // InputMHFEM
template void cuDFNsys::InputObjectData<double>::InputMHFEM(const string &filename,
                                                            cuDFNsys::MHFEM<double> &MHFEM);
template void cuDFNsys::InputObjectData<float>::InputMHFEM(const string &filename,
                                                           cuDFNsys::MHFEM<float> &MHFEM);