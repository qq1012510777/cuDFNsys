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