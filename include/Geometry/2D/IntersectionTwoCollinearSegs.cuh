///////////////////////////////////////////////////////////////////
// NAME:              IntersectionTwoCollinearSegs.cuh
//
// PURPOSE:           Identify intersection between two collinear segments
//
// FUNCTIONS/OBJECTS: IntersectionTwoCollinearSegs
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#include "../../GlobalDef/GlobalDef.cuh"

namespace cuDFNsys
{
__device__ __host__ bool IntersectionTwoCollinearSegs(float2 *Seg_1,
                                                      float2 *Seg_2,
                                                      float2 *intersection,
                                                      int *sign,
                                                      float _TOL_);
}; // namespace cuDFNsys