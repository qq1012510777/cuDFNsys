///////////////////////////////////////////////////////////////////
// NAME:              IfParticlePositionNearOneVertexOfElement.cuh
//
// PURPOSE:           If the particle position is very close to one of the vertexes of grid
//
// FUNCTIONS/OBJECTS: IfParticlePositionNearOneVertexOfElement
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../DataTypeSelector/DataTypeSelector.cuh"
#include "../GlobalDef/GlobalDef.cuh"

namespace cuDFNsys
{
template <typename T>
__device__ __host__ bool IfParticlePositionNearOneVertexOfElement(cuDFNsys::Vector2<T> Position_p, cuDFNsys::Vector2<T> Tri[3], T _TOL_p)
{
    for (uint i = 0; i < 3; ++i)
    {
        cuDFNsys::Vector2<T> V;
        V.x = Position_p.x - Tri[i].x;
        V.y = Position_p.y - Tri[i].y;

        T norm_ = sqrt(V.x * V.x + V.y * V.y);

        if (norm_ < _TOL_p)
            return true;
    }
    return false;
};
template __device__ __host__ bool IfParticlePositionNearOneVertexOfElement<double>(cuDFNsys::Vector2<double> Position_p, cuDFNsys::Vector2<double> Tri[3], double _TOL_p);
template __device__ __host__ bool IfParticlePositionNearOneVertexOfElement<float>(cuDFNsys::Vector2<float> Position_p, cuDFNsys::Vector2<float> Tri[3], float _TOL_p);
}; // namespace cuDFNsys