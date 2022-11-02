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
// NAME:              EdgeToEle.cuh
//
// PURPOSE:           a struct of EdgeToEle: record the shared element IDs of each
//                    (separated) edge
//
// FUNCTIONS/OBJECTS: EdgeToEle
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../GlobalDef/GlobalDef.cuh"
#include "../DataTypeSelector/DataTypeSelector.cuh"

namespace cuDFNsys
{
struct EdgeToEle
{
    // number of shared elements for an edge, (at least = 1, at most = _NumOfSharedEleAtMost in GlobalDef.cuh)
    uint NumSharedEle = 0;
    // here, the default of number of shared elements is _NumOfSharedEleAtMost.
    uint EleID[_NumOfSharedEleAtMost] = {0};
    // the shared edge (actually the same edge, but different local NO. values); local edge NO: 0, 1, or 2
    uint LocalEdgeNO[_NumOfSharedEleAtMost] = {0};
};
}; // namespace cuDFNsys