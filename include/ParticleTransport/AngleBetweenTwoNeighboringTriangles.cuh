///////////////////////////////////////////////////////////////////
// NAME:              AngleBetweenTwoNeighboringTriangles.cuh
//
// PURPOSE:           return the angle between two neighboring triangles
//
// FUNCTIONS/OBJECTS: AngleBetweenTwoNeighboringTriangles
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../DataTypeSelector/DataTypeSelector.cuh"
#include "../GlobalDef/GlobalDef.cuh"
#include "../MatrixManipulation/MatrixManipulation.cuh"

namespace cuDFNsys
{
template <typename T>
__host__ __device__ T AngleBetweenTwoNeighboringTriangles(const cuDFNsys::Vector3<T> preEle[3],
                                                          const cuDFNsys::Vector3<T> nextEle[3],
                                                          const uint &localEdgeNo_preEle,
                                                          const uint &localEdgeNo_nextEle/*,
                                                          cuDFNsys::Vector3<T> &d1,
                                                          cuDFNsys::Vector3<T> &d2*/)
{
    cuDFNsys::Vector3<T> pntA = preEle[(localEdgeNo_preEle + 2) % 3];
    cuDFNsys::Vector3<T> pntB = preEle[localEdgeNo_preEle];
    cuDFNsys::Vector3<T> pntC = preEle[(localEdgeNo_preEle + 1) % 3];

    cuDFNsys::Vector3<T> pntD = nextEle[(localEdgeNo_nextEle + 2) % 3];

    // printf("pntA: %.40f, %.40f,  %.40f\n", pntA.x, pntA.y, pntA.z);
    // printf("pntB: %.40f, %.40f,  %.40f\n", pntB.x, pntB.y, pntB.z);
    // printf("pntC: %.40f, %.40f,  %.40f\n", pntC.x, pntC.y, pntC.z);
    // printf("pntD: %.40f, %.40f,  %.40f\n", pntD.x, pntD.y, pntD.z);

    cuDFNsys::Vector3<T> V = cuDFNsys::MakeVector3(pntC.x - pntB.x, pntC.y - pntB.y, pntC.z - pntB.z);
    T norm_V = sqrt(V.x * V.x + V.y * V.y + V.z * V.z);
    V.x /= norm_V;
    V.y /= norm_V;
    V.z /= norm_V;

    cuDFNsys::Vector3<T> BA = cuDFNsys::MakeVector3(pntA.x - pntB.x, pntA.y - pntB.y, pntA.z - pntB.z);
    T norm_BA = sqrt(BA.x * BA.x + BA.y * BA.y + BA.z * BA.z);
    BA.x /= norm_BA;
    BA.y /= norm_BA;
    BA.z /= norm_BA;

    cuDFNsys::Vector3<T> BD = cuDFNsys::MakeVector3(pntD.x - pntB.x, pntD.y - pntB.y, pntD.z - pntB.z);
    T norm_BD = sqrt(BD.x * BD.x + BD.y * BD.y + BD.z * BD.z);
    BD.x /= norm_BD;
    BD.y /= norm_BD;
    BD.z /= norm_BD;

    cuDFNsys::Vector3<T> a = cuDFNsys::CrossProductVector3<T>(BA, V);
    a = cuDFNsys::CrossProductVector3<T>(V, a);
    T norm_a = sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
    a.x /= norm_a;
    a.y /= norm_a;
    a.z /= norm_a;
    //d1 = a;

    cuDFNsys::Vector3<T> b = cuDFNsys::CrossProductVector3<T>(BD, V);
    b = cuDFNsys::CrossProductVector3<T>(V, b);
    T norm_b = sqrt(b.x * b.x + b.y * b.y + b.z * b.z);
    b.x /= norm_b;
    b.y /= norm_b;
    b.z /= norm_b;
    //d2 = b;

    double theta = 0;
    theta = acos((a.x * b.x + a.y * b.y + a.z * b.z));
    return theta;
};
template __host__ __device__ double AngleBetweenTwoNeighboringTriangles<double>(const cuDFNsys::Vector3<double> preEle[3],
                                                                                const cuDFNsys::Vector3<double> nextEle[3],
                                                                                const uint &localEdgeNo_preEle,
                                                                                const uint &localEdgeNo_nextEle/*,
                                                                                cuDFNsys::Vector3<double> &d1,
                                                                                cuDFNsys::Vector3<double> &d2*/);
template __host__ __device__ float AngleBetweenTwoNeighboringTriangles<float>(const cuDFNsys::Vector3<float> preEle[3],
                                                                              const cuDFNsys::Vector3<float> nextEle[3],
                                                                              const uint &localEdgeNo_preEle,
                                                                              const uint &localEdgeNo_nextEle/*,
                                                                              cuDFNsys::Vector3<float> &d1,
                                                                              cuDFNsys::Vector3<float> &d2*/);
}; // namespace cuDFNsys