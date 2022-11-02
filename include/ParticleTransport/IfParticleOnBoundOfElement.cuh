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
// NAME:              IfParticleOnBoundOfElement.cuh
//
// PURPOSE:           If particle lies on bound of element
//
// FUNCTIONS/OBJECTS: IfParticleOnBoundOfElement
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../DataTypeSelector/DataTypeSelector.cuh"
#include "../Geometry/Geometry.cuh"
#include "../GlobalDef/GlobalDef.cuh"

namespace cuDFNsys
{
template <typename T>
__host__ __device__ bool IfParticleOnBoundOfElement(cuDFNsys::Vector2<T> PositionP,
                                                    cuDFNsys::Vector2<T> GridVertex[3],
                                                    int &EdgeNOLocal,
                                                    T _TOL_)
{

    for (uint i = 0; i < 3; ++i)
    {
        int o1 = cuDFNsys::OrientationThree2DPnts<T>(GridVertex[i],
                                                     GridVertex[(i + 1) % 3],
                                                     PositionP, _TOL_);
        // printf("In 'IfParticleOnBoundOfElement', cuDFNsys::OrientationThree2DPnts: %d\n", o1);
        if (o1 == 0 && cuDFNsys::If2DPntLiesOnCollinearSeg<T>(GridVertex[i],
                                                              PositionP,
                                                              GridVertex[(i + 1) % 3]))
        {
            EdgeNOLocal = (int)i;
            return true;
        };
    };

    EdgeNOLocal = -1;
    return false;
};
template __host__ __device__ bool IfParticleOnBoundOfElement<double>(cuDFNsys::Vector2<double> PositionP,
                                                                     cuDFNsys::Vector2<double> GridVertex[3],
                                                                     int &EdgeNOLocal,
                                                                     double _TOL_);
template __host__ __device__ bool IfParticleOnBoundOfElement<float>(cuDFNsys::Vector2<float> PositionP,
                                                                    cuDFNsys::Vector2<float> GridVertex[3],
                                                                    int &EdgeNOLocal,
                                                                    float _TOL_);
}; // namespace cuDFNsys