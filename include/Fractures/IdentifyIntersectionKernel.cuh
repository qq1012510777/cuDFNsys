///////////////////////////////////////////////////////////////////
// NAME:              IdentifyIntersectionKernel.cuh
//
// PURPOSE:           Identify intersection using GPU
//
// FUNCTIONS/OBJECTS: IdentifyIntersectionKernel
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../Geometry/2D/Intersection2DLine2DPoly.cuh"
#include "../Geometry/2D/IntersectionTwoCollinearSegs.cuh"
#include "../Geometry/3D/Intersection3DPolyXYPlane.cuh"
#include "../GlobalDef/GlobalDef.cuh"
#include "../MatrixManipulation/MatrixManipulation.cuh"
#include "Fracture.cuh"
#include "Intersection.cuh"
using namespace std;

namespace cuDFNsys
{
__global__ void IdentifyIntersectionKernel(cuDFNsys::Fracture *verts,
                                           int count,
                                           cuDFNsys::Intersection *Int_sec,
                                           bool if_trucncated);
}; // namespace cuDFNsys