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
// NAME:              Geometry.cuh
//
// PURPOSE:           The API of Geometry
//
// FUNCTIONS/OBJECTS: N/A
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////

#pragma once
#include "../DataTypeSelector/DataTypeSelector.cuh"
#include "../GlobalDef/GlobalDef.cuh"
#include "./2D/DistancePnt2DSeg.cuh"
#include "./2D/If2DPntLiesOnCollinearSeg.cuh"
#include "./2D/IfPntInside2DConvexPoly.cuh"
#include "./2D/IfPntLiesOnBound2DConvexPoly.cuh"
#include "./2D/Intersection2DLine2DPoly.cuh"
#include "./2D/Intersection2DLine2DSeg.cuh"
#include "./2D/IntersectionTwoCollinear2DSegs.cuh"
#include "./2D/OrientationThree2DPnts.cuh"
#include "./2D/Scale2DSegment.cuh"
#include "./2D/Triangle2DArea.cuh"
#include "./2D/Triangle2DOrientation.cuh"
#include "./3D/DistancePnt3DPlane.cuh"
#include "./3D/If3DTriangleSkinny.cuh"
#include "./3D/Intersection3DPolyXYPlane.cuh"
#include "./3D/Intersection3DSegXYPlane.cuh"
#include "./3D/Scale3DTriangle.cuh"
#include "./3D/Triangle3DArea.cuh"