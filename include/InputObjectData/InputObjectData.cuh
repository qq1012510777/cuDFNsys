///////////////////////////////////////////////////////////////////
// NAME:              InputObjectData.cuh
//
// PURPOSE:           Input data of fractures, mesh, and so on
//
// FUNCTIONS/OBJECTS: InputObjectData
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
class InputObjectData
{
public:
    // constructor
    InputObjectData(){};

public:
    void InputFractures(const string &filename,
                        thrust::host_vector<cuDFNsys::Fracture<T>> &Frac_verts_host);

    // input mesh data
    void InputMesh(const string &filename,
                   cuDFNsys::Mesh<T> &mesh,
                   std::vector<size_t> *Fracs_percol);
};
}; // namespace cuDFNsys