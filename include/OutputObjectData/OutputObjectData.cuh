///////////////////////////////////////////////////////////////////
// NAME:              OutputObjectData.cuh
//
// PURPOSE:           Output data of fractures, mesh, and so on
//
// FUNCTIONS/OBJECTS: OutputObjectData
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../DataTypeSelector/DataTypeSelector.cuh"
#include "../Fractures/Fracture.cuh"
#include "../GlobalDef/GlobalDef.cuh"
#include "../HDF5API/HDF5API.cuh"
#include "../Mesh/Mesh.cuh"

namespace cuDFNsys
{
template <typename T>
class OutputObjectData
{
public:
    // constructor
    OutputObjectData(){};

public:
    // output Fractures in H5
    void OutputFractures(const string &filename_, const thrust::host_vector<cuDFNsys::Fracture<T>> &Frac_verts_host);

    // output Mesh in h5
    void OutputMesh(const string &filename_,
                    cuDFNsys::Mesh<T> mesh,
                    const std::vector<size_t> &Fracs_percol);
};
}; // namespace cuDFNsys