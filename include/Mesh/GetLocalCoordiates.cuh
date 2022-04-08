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
#include "../Fractures/Fracture.cuh"
#include "../Geometry/2D/Triangle2DOrientation.cuh"
#include "../GlobalDef/GlobalDef.cuh"
#include "../MatrixManipulation/MatrixManipulation.cuh"
#include "EleCoor.cuh"

namespace cuDFNsys
{
__global__ void GetLocalCoordiates(uint3 *element_3D_dev_ptr,
                                   cuDFNsys::Fracture *Frac_verts_device_ptr,
                                   uint *element_Frac_Tag_dev_ptr,
                                   cuDFNsys::EleCoor *coordinate_2D_dev_ptr,
                                   float3 *coordinate_3D_dev_ptr,
                                   int ele_count);
}; // namespace cuDFNsys