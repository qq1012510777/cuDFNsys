///////////////////////////////////////////////////////////////////
// NAME:              GetLocalCoordiates.cuh
//
// PURPOSE:           Get 2D local coordiates of elements
//
// FUNCTIONS/OBJECTS: GetLocalCoordiates
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../DataTypeSelector/DataTypeSelector.cuh"
#include "../Fractures/Fracture.cuh"
#include "../Geometry/Geometry.cuh"
#include "../GlobalDef/GlobalDef.cuh"
#include "../MatrixManipulation/MatrixManipulation.cuh"
#include "EleCoor.cuh"

namespace cuDFNsys
{
template <typename T>
__global__ void GetLocalCoordiates(uint3 *element_3D_dev_ptr,
                                   cuDFNsys::Fracture<T> *Frac_verts_device_ptr,
                                   uint *element_Frac_Tag_dev_ptr,
                                   cuDFNsys::EleCoor<T> *coordinate_2D_dev_ptr,
                                   cuDFNsys::Vector3<T> *coordinate_3D_dev_ptr,
                                   int ele_count);
}; // namespace cuDFNsys