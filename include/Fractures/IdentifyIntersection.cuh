///////////////////////////////////////////////////////////////////
// NAME:              IdentifyIntersection.cuh
//
// PURPOSE:           Identify intersection
//
// FUNCTIONS/OBJECTS: IdentifyIntersection
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
#include "IdentifyFracPairSphericalDetection.cuh"
#include "IdentifyIntersectionKernel.cuh"
#include "Intersection.cuh"
#include <unistd.h>

namespace cuDFNsys
{
class IdentifyIntersection
{
public:
    // constructor CPU
    IdentifyIntersection(thrust::host_vector<cuDFNsys::Fracture> verts,
                         const bool &if_trucncated,
                         MapIntersection &Intersection_map);
    // constructor GPU
    IdentifyIntersection(const size_t &Fracsize,
                         cuDFNsys::Fracture *Frac_verts_device_ptr,
                         const bool &if_trucncated,
                         MapIntersection &Intersection_map);
};
}; // namespace cuDFNsys