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