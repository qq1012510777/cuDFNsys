///////////////////////////////////////////////////////////////////
// NAME:              GlobalDef.cuh
//
// PURPOSE:           Define some global variables
//
// FUNCTIONS/OBJECTS: N/A
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////

#pragma once
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <cuda.h>
#include <cuda_runtime.h>
#include <curand.h>
#include <curand_kernel.h>
#include <device_launch_parameters.h>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <sstream>
#include <stdlib.h>
#include <sys/time.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <vector>

using namespace std;

typedef std::map<pair<size_t, size_t>, pair<float3, float3>> MapIntersection;

#define _TOL_Intersection3DPolyXYPlane 1e-7
#define _TOL_Intersection2DLine2DPoly 1e-7
#define _TOL_IntersectionTwoCollinearSegs 1e-7
#define _TOL_If3DTriangleSkinny 1e-5
#define _TOL_IdentifyEleFrac 1e-2
#define _TOL_IfTwoEndsDirchlet 1e-7

__device__ const float _TOL_ParticleOnGridBound = 1e-7;

const uint _NumOfSharedEleAtMost = 4;