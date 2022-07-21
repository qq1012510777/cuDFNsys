///////////////////////////////////////////////////////////////////
// NAME:              Roate2DPositionTo3D.cuh
//
// PURPOSE:           Rotate a 2D particle position to 3D
//
// FUNCTIONS/OBJECTS: Roate2DPositionTo3D
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../DataTypeSelector/DataTypeSelector.cuh"
#include "../Fractures/Fracture.cuh"
#include "../GlobalDef/GlobalDef.cuh"
#include "../MatrixManipulation/MatrixManipulation.cuh"

namespace cuDFNsys
{
template <typename T>
__host__ __device__ cuDFNsys::Vector3<T> Roate2DPositionTo3D(cuDFNsys::Vector2<T> PositionP,
                                                             cuDFNsys::Fracture<T> OneFrac)
{
    cuDFNsys::Vector3<T> P_3D = cuDFNsys::MakeVector3(PositionP.x, PositionP.y, (T)0.0);

    T RK_2[3][3];
    OneFrac.RoationMatrix(RK_2, 23);

    P_3D = cuDFNsys::ProductSquare3Vector3<T>(RK_2, P_3D);
    P_3D.x += OneFrac.Center.x;
    P_3D.y += OneFrac.Center.y;
    P_3D.z += OneFrac.Center.z;

    return P_3D;
};
template __host__ __device__ cuDFNsys::Vector3<double> Roate2DPositionTo3D<double>(cuDFNsys::Vector2<double> PositionP,
                                                                                   cuDFNsys::Fracture<double> OneFrac);
template __host__ __device__ cuDFNsys::Vector3<float> Roate2DPositionTo3D<float>(cuDFNsys::Vector2<float> PositionP,
                                                                                 cuDFNsys::Fracture<float> OneFrac);
}; // namespace cuDFNsys