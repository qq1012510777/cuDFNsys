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

#include "OutputObjectData/OutputObjectData.cuh"

// ====================================================
// NAME:        OutputFractures
// DESCRIPTION: OutputFractures h5 file
// AUTHOR:      Tingchang YIN
// DATE:        02/08/2022
// ====================================================
template <typename T>
void cuDFNsys::OutputObjectData<T>::OutputFractures(
    const string &filename_,
    const thrust::host_vector<cuDFNsys::Fracture<T>> &Frac_verts_host,
    const T &L, double3 DomainDimensionRatio)
{
    uint Dsize = Frac_verts_host.size();

    vector<string> datasername_ = {
        "Conductivity",      "Verts3D", "Center",           "Verts3DTruncated",
        "NumVertsTruncated", "Radius",  "ConnectModelSurf", "NormalVec"};

    cuDFNsys::HDF5API h5_;
    h5_.NewFile(filename_);

    uint2 dim_e = make_uint2(1, 1);
    uint DiszeP[1] = {Dsize};
    h5_.AddDataset(filename_, "N", "NumFractures", DiszeP, dim_e);

    uint2 dim_q = make_uint2(3, 1);
    double sd[3] = {DomainDimensionRatio.x, DomainDimensionRatio.y,
                    DomainDimensionRatio.z};
    h5_.AddDataset(filename_, "N", "DomainDimensionRatio", sd, dim_q);

    T L_m[1] = {L};
    h5_.AddDataset(filename_, "N", "L", L_m, dim_e);

    for (uint i = 0; i < Dsize; ++i)
    {
        string groupname_ = "Fracture_" + std::to_string(i + 1);

        T Conductivity_[1] = {Frac_verts_host[i].Conductivity};
        T Verts3D_[12] = {
            Frac_verts_host[i].Verts3D[0].x, Frac_verts_host[i].Verts3D[1].x,
            Frac_verts_host[i].Verts3D[2].x, Frac_verts_host[i].Verts3D[3].x,
            Frac_verts_host[i].Verts3D[0].y, Frac_verts_host[i].Verts3D[1].y,
            Frac_verts_host[i].Verts3D[2].y, Frac_verts_host[i].Verts3D[3].y,
            Frac_verts_host[i].Verts3D[0].z, Frac_verts_host[i].Verts3D[1].z,
            Frac_verts_host[i].Verts3D[2].z, Frac_verts_host[i].Verts3D[3].z};

        T Center_[3] = {Frac_verts_host[i].Center.x,
                        Frac_verts_host[i].Center.y,
                        Frac_verts_host[i].Center.z};
        T Verts3DTruncated_[24] = {Frac_verts_host[i].Verts3DTruncated[0].x,
                                   Frac_verts_host[i].Verts3DTruncated[1].x,
                                   Frac_verts_host[i].Verts3DTruncated[2].x,
                                   Frac_verts_host[i].Verts3DTruncated[3].x,
                                   Frac_verts_host[i].Verts3DTruncated[4].x,
                                   Frac_verts_host[i].Verts3DTruncated[5].x,
                                   Frac_verts_host[i].Verts3DTruncated[6].x,
                                   Frac_verts_host[i].Verts3DTruncated[7].x,
                                   Frac_verts_host[i].Verts3DTruncated[0].y,
                                   Frac_verts_host[i].Verts3DTruncated[1].y,
                                   Frac_verts_host[i].Verts3DTruncated[2].y,
                                   Frac_verts_host[i].Verts3DTruncated[3].y,
                                   Frac_verts_host[i].Verts3DTruncated[4].y,
                                   Frac_verts_host[i].Verts3DTruncated[5].y,
                                   Frac_verts_host[i].Verts3DTruncated[6].y,
                                   Frac_verts_host[i].Verts3DTruncated[7].y,
                                   Frac_verts_host[i].Verts3DTruncated[0].z,
                                   Frac_verts_host[i].Verts3DTruncated[1].z,
                                   Frac_verts_host[i].Verts3DTruncated[2].z,
                                   Frac_verts_host[i].Verts3DTruncated[3].z,
                                   Frac_verts_host[i].Verts3DTruncated[4].z,
                                   Frac_verts_host[i].Verts3DTruncated[5].z,
                                   Frac_verts_host[i].Verts3DTruncated[6].z,
                                   Frac_verts_host[i].Verts3DTruncated[7].z};

        T NumVertsTruncated_[1] = {(T)Frac_verts_host[i].NumVertsTruncated};
        T Radius_[1] = {Frac_verts_host[i].Radius};
        T ConnectModelSurf_[6] = {(T)Frac_verts_host[i].ConnectModelSurf[0],
                                  (T)Frac_verts_host[i].ConnectModelSurf[1],
                                  (T)Frac_verts_host[i].ConnectModelSurf[2],
                                  (T)Frac_verts_host[i].ConnectModelSurf[3],
                                  (T)Frac_verts_host[i].ConnectModelSurf[4],
                                  (T)Frac_verts_host[i].ConnectModelSurf[5]};
        T NormalVec_[3] = {Frac_verts_host[i].NormalVec.x,
                           Frac_verts_host[i].NormalVec.y,
                           Frac_verts_host[i].NormalVec.z};

        // vector<string> datasername_ = {
        //     "Conductivity",
        //     "Verts3D",
        //     "Center",
        //     "Verts3DTruncated",
        //     "NumVertsTruncated",
        //     "Radius",
        //     "ConnectModelSurf",
        //     "NormalVec"};

        vector<T *> data = {Conductivity_,     Verts3D_,           Center_,
                            Verts3DTruncated_, NumVertsTruncated_, Radius_,
                            ConnectModelSurf_, NormalVec_};
        vector<uint2> dim = {make_uint2(1, 1), make_uint2(3, 4),
                             make_uint2(1, 3), make_uint2(3, 8),
                             make_uint2(1, 1), make_uint2(1, 1),
                             make_uint2(1, 6), make_uint2(1, 3)};
        h5_.AddDatasetsWithOneGroup(filename_, groupname_, datasername_, data,
                                    dim);
    }
}; // OutputFractures
template void cuDFNsys::OutputObjectData<double>::OutputFractures(
    const string &filename_,
    const thrust::host_vector<cuDFNsys::Fracture<double>> &Frac_verts_host,
    const double &L, double3 DomainDimensionRatio);
template void cuDFNsys::OutputObjectData<float>::OutputFractures(
    const string &filename_,
    const thrust::host_vector<cuDFNsys::Fracture<float>> &Frac_verts_host,
    const float &L, double3 DomainDimensionRatio);

// ====================================================
// NAME:        OutputMesh
// DESCRIPTION: OutputMesh h5 file
// AUTHOR:      Tingchang YIN
// DATE:        02/08/2022
// ====================================================
template <typename T>
void cuDFNsys::OutputObjectData<T>::OutputMesh(
    const string &filename_, cuDFNsys::Mesh<T> mesh,
    const std::vector<size_t> &Fracs_percol)
{
    cuDFNsys::HDF5API h5_;
    h5_.NewFile(filename_);

    uint2 dim_e = make_uint2(1, 1);
    uint NumFractures[1] = {(uint)Fracs_percol.size()};
    h5_.AddDataset(filename_, "N", "NumFractures", NumFractures, dim_e);

    h5_.AddDataset(filename_, "N", "MeanGridSize", &mesh.MeanGridSize, dim_e);

    uint NumNodes[1] = {(uint)mesh.Coordinate3D.size()};
    h5_.AddDataset(filename_, "N", "NumNodes", NumNodes, dim_e);

    uint NumElements[1] = {(uint)mesh.Element3D.size()};
    h5_.AddDataset(filename_, "N", "NumElements", NumElements, dim_e);

    uint NumInletEdge[1] = {(uint)mesh.InletEdgeNOLen.size()};
    h5_.AddDataset(filename_, "N", "NumInletEdges", NumInletEdge, dim_e);

    uint NumOutletEdge[1] = {(uint)mesh.OutletEdgeNOLen.size()};
    h5_.AddDataset(filename_, "N", "NumOutletEdges", NumOutletEdge, dim_e);

    std::vector<uint> Fracs_percol_II(Fracs_percol.size());
    std::copy(Fracs_percol.begin(), Fracs_percol.end(), Fracs_percol_II.data());

    h5_.AddDataset<uint>(filename_, "N", "Fracs_percol", Fracs_percol_II.data(),
                         make_uint2(Fracs_percol.size(), 0));

    h5_.AddDataset<T>(filename_, "N", "InletTraceLength",
                      &(mesh.InletTraceLength), make_uint2(1, 1));
    h5_.AddDataset<T>(filename_, "N", "OutletTraceLength",
                      &(mesh.OutletTraceLength), make_uint2(1, 1));

    vector<string> datasetname = {"MeshSuccess",     "FracID",
                                  "Coordinate3D",    "Element2D",
                                  "NumElementFrac",  "Element3D",
                                  "Coordinate2D",    "ElementFracTag",
                                  "EdgeAttri",       "InletEdgeNOLen",
                                  "OutletEdgeNOLen", "NumInteriorEdges",
                                  "NumNeumannEdges", "Dir"};

    T MeshSuccess[1] = {(T)mesh.MeshSuccess};

    T *FracID = new T[NumFractures[0]];
    for (uint i = 0; i < NumFractures[0]; ++i)
        FracID[i] = (T)Fracs_percol[i];

    T *Coordinate3D = new T[3 * NumNodes[0]];
    for (uint i = 0; i < 3; ++i)
        for (uint j = 0; j < NumNodes[0]; ++j)
        {
            T *KL = &(mesh.Coordinate3D[j].x);
            Coordinate3D[j + i * NumNodes[0]] = KL[i];
        }

    T *Element2D = new T[3 * NumElements[0]];
    T *NumElementFrac = new T[NumFractures[0]];

    uint iio = 0;
    for (uint i = 0; i < 3; ++i)
        for (uint j = 0; j < NumFractures[0]; ++j)
            for (uint k = 0; k < mesh.Element2D[j].size(); ++k)
            {
                NumElementFrac[j] = mesh.Element2D[j].size();

                uint *KL = &(mesh.Element2D[j][k].x);
                Element2D[iio] = (T)KL[i];
                iio++;
            };

    T *Element3D = new T[3 * NumElements[0]];
    for (uint i = 0; i < 3; ++i)
        for (uint j = 0; j < NumElements[0]; ++j)
        {
            uint *KL = &(mesh.Element3D[j].x);
            Element3D[i * NumElements[0] + j] = (T)KL[i];
        }

    T *Coordinate2D = new T[NumElements[0] * 6];
    for (uint i = 0; i < 6; ++i)
        for (uint j = 0; j < NumElements[0]; ++j)
        {
            T *KL;

            if (i < 3)
                KL = &(mesh.Coordinate2D[j].x[0]);
            else
                KL = &(mesh.Coordinate2D[j].y[0]);

            Coordinate2D[i * NumElements[0] + j] = KL[i % 3];
        }

    T *ElementFracTag = new T[NumElements[0]];
    for (uint i = 0; i < NumElements[0]; ++i)
        ElementFracTag[i] = (T)mesh.ElementFracTag[i];

    T *EdgeAttri = new T[6 * NumElements[0]];
    for (uint i = 0; i < 6; ++i)
        for (uint j = 0; j < NumElements[0]; ++j)
        {
            int *KL;

            if (i < 3)
                KL = &(mesh.EdgeAttri[j].e[0]);
            else
                KL = &(mesh.EdgeAttri[j].no[0]);

            EdgeAttri[i * NumElements[0] + j] = (T)KL[i % 3];
        }

    T *InletEdgeNOLen = new T[NumInletEdge[0] * 2];
    for (uint i = 0; i < 2; ++i)
        for (uint j = 0; j < NumInletEdge[0]; ++j)
        {
            T *KL = &(mesh.InletEdgeNOLen[j].x);
            InletEdgeNOLen[i * NumInletEdge[0] + j] = KL[i];
        }

    T *OutletEdgeNOLen = new T[NumOutletEdge[0] * 2];
    for (uint i = 0; i < 2; ++i)
        for (uint j = 0; j < NumOutletEdge[0]; ++j)
        {
            T *KL = &(mesh.OutletEdgeNOLen[j].x);
            OutletEdgeNOLen[i * NumOutletEdge[0] + j] = KL[i];
        }

    T NumInteriorEdges[1] = {(T)mesh.NumInteriorEdges};

    T NumNeumannEdges[1] = {(T)mesh.NumNeumannEdges};

    T Dir[1] = {(T)mesh.Dir};

    vector<uint2> dim = {make_uint2(1, 1),
                         make_uint2(NumFractures[0], 1),
                         make_uint2(3, NumNodes[0]),
                         make_uint2(3, NumElements[0]),
                         make_uint2(NumFractures[0], 1),
                         make_uint2(3, NumElements[0]),
                         make_uint2(6, NumElements[0]),
                         make_uint2(NumElements[0], 1),
                         make_uint2(6, NumElements[0]),
                         make_uint2(2, NumInletEdge[0]),
                         make_uint2(2, NumOutletEdge[0]),
                         make_uint2(1, 1),
                         make_uint2(1, 1),
                         make_uint2(1, 1)};

    vector<T *> data = {MeshSuccess,     FracID,
                        Coordinate3D,    Element2D,
                        NumElementFrac,  Element3D,
                        Coordinate2D,    ElementFracTag,
                        EdgeAttri,       InletEdgeNOLen,
                        OutletEdgeNOLen, NumInteriorEdges,
                        NumNeumannEdges, Dir};

    h5_.AddDatasetsWithOneGroup(filename_, "group_mesh", datasetname, data,
                                dim);

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
}; // OutputMesh
template void cuDFNsys::OutputObjectData<double>::OutputMesh(
    const string &filename_, cuDFNsys::Mesh<double> mesh,
    const std::vector<size_t> &Fracs_percol);
template void cuDFNsys::OutputObjectData<float>::OutputMesh(
    const string &filename_, cuDFNsys::Mesh<float> mesh,
    const std::vector<size_t> &Fracs_percol);

// ====================================================
// NAME:        OutputMHFEM
// DESCRIPTION: OutputMHFEM h5 file
// AUTHOR:      Tingchang YIN
// DATE:        27/12/2022
// ====================================================
template <typename T>
void cuDFNsys::OutputObjectData<T>::OutputMHFEM(const string &filename_,
                                                cuDFNsys::MHFEM<T> mhfem)
{
    //
    cuDFNsys::HDF5API h5_;
    h5_.NewFile(filename_);

    uint2 dim_e = make_uint2(1, 1);
    T Value[1] = {mhfem.QIn};
    h5_.AddDataset(filename_, "N", "QIn", Value, dim_e);

    Value[0] = mhfem.QOut;
    h5_.AddDataset(filename_, "N", "QOut", Value, dim_e);

    Value[0] = mhfem.MuOverRhoG;
    h5_.AddDataset(filename_, "N", "MuOverRhoG", Value, dim_e);

    Value[0] = (mhfem.IfPeriodic ? 1 : 0);
    h5_.AddDataset(filename_, "N", "IfPeriodic", Value, dim_e);

    Value[0] = mhfem.ConsTq;
    h5_.AddDataset(filename_, "N", "ConsTq", Value, dim_e);

    Value[0] = mhfem.InletLength;
    h5_.AddDataset(filename_, "N", "InletLength", Value, dim_e);

    Value[0] = mhfem.OutletLength;
    h5_.AddDataset(filename_, "N", "OutletLength", Value, dim_e);

    Value[0] = mhfem.QError;
    h5_.AddDataset(filename_, "N", "QError", Value, dim_e);

    Value[0] = mhfem.Permeability;
    h5_.AddDataset(filename_, "N", "Permeability", Value, dim_e);

    Value[0] = mhfem.InletP;
    h5_.AddDataset(filename_, "N", "InletP", Value, dim_e);

    Value[0] = mhfem.OutletP;
    h5_.AddDataset(filename_, "N", "OutletP", Value, dim_e);

    double Value2[1] = {mhfem.TripletTime};
    h5_.AddDataset(filename_, "N", "TripletTime", Value2, dim_e);

    h5_.AddDataset(filename_, "N", "MeanVelocity", &mhfem.MeanVelocity, dim_e);
    h5_.AddDataset(filename_, "N", "MaxVelocity", &mhfem.MaxVelocity, dim_e);

    int Value3[1] = {mhfem.Dir};
    h5_.AddDataset(filename_, "N", "Dir", Value3, dim_e);

    vector<double> PressureInteriorEdge_(mhfem.PressureInteriorEdge.data(),
                                         mhfem.PressureInteriorEdge.data() +
                                             mhfem.PressureInteriorEdge.rows());
    uint2 dim_r = make_uint2(1, PressureInteriorEdge_.size());
    h5_.AddDataset(filename_, "N", "PressureInteriorEdge",
                   PressureInteriorEdge_.data(), dim_r);

    vector<double> PressureEles_(mhfem.PressureEles.data(),
                                 mhfem.PressureEles.data() +
                                     mhfem.PressureEles.rows());
    dim_r = make_uint2(1, PressureEles_.size());
    h5_.AddDataset(filename_, "N", "PressureEles", PressureEles_.data(), dim_r);

    vector<double> VelocityNormalScalarSepEdges_(
        mhfem.VelocityNormalScalarSepEdges.data(),
        mhfem.VelocityNormalScalarSepEdges.data() +
            mhfem.VelocityNormalScalarSepEdges.rows());
    dim_r = make_uint2(1, VelocityNormalScalarSepEdges_.size());
    h5_.AddDataset(filename_, "N", "VelocityNormalScalarSepEdges",
                   VelocityNormalScalarSepEdges_.data(), dim_r);

}; // OutputMHFEM
template void
cuDFNsys::OutputObjectData<double>::OutputMHFEM(const string &filename_,
                                                cuDFNsys::MHFEM<double> mhfem);
template void
cuDFNsys::OutputObjectData<float>::OutputMHFEM(const string &filename_,
                                               cuDFNsys::MHFEM<float> mhfem);