///////////////////////////////////////////////////////////////////
// NAME:              IdentifyEleFrac.cuh
//
// PURPOSE:           Identify fracture ID of elements
//
// FUNCTIONS/OBJECTS: IdentifyEleFrac
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../Fractures/Fracture.cuh"
#include "../Geometry/2D/IfPntInside2DConvexPoly.cuh"
#include "../Geometry/2D/IfPntLiesOnBound2DConvexPoly.cuh"
#include "../Geometry/3D/DistancePnt3DPlane.cuh"
#include "../GlobalDef/GlobalDef.cuh"

namespace cuDFNsys
{
__global__ void IdentifyEleFrac(uint3 *One_entity_one_ele_dev_ptr,
                                float3 *coordinate_3D_dev_ptr,
                                cuDFNsys::Fracture *Frac_verts_device_ptr,
                                int *Elements_Frac_dev_ptr,
                                int entity_count,
                                int frac_count,
                                float _tol_);
}; // namespace cuDFNsys