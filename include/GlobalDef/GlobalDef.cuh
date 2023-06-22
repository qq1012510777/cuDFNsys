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
#include <thrust/sequence.h>
#include <vector>

using namespace std;

#define _TOL_Intersection3DPolyXYPlane 1e-7
#define _TOL_Intersection2DLine2DPoly 1e-7
#define _TOL_IntersectionTwoCollinearSegs 1e-7
#define _TOL_If3DTriangleSkinny 1e-5
#define _TOL_IdentifyEleFrac 1e-3
#define _TOL_IfTwoEndsDirchlet 1e-7

__device__ const float _TOL_ParticleOnGridBound = 1e-3;
__device__ const uint _SizeOfArray_CrossedGlobalEdge_ = 50;

const uint _NumOfSharedEleAtMost = 4;
const uint _NumOfNeighborEleAtMost = 80;
const uint _ParTran_MaxLoopTimes = 50;