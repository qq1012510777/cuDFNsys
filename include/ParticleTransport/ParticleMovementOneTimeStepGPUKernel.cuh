///////////////////////////////////////////////////////////////////
// NAME:              ParticleMovementOneTimeStepGPUKernel.cuh
//
// PURPOSE:           GPU kernel function to move particles
//                    in one time step
//
// FUNCTIONS/OBJECTS: ParticleMovementOneTimeStepGPUKernel
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "../Geometry/2D/DistancePnt2DSeg.cuh"
#include "../Geometry/2D/IfPntInside2DConvexPoly.cuh"
#include "../Geometry/2D/IfPntLiesOnBound2DConvexPoly.cuh"
#include "../Geometry/2D/IntersectionTwo2DSegs.cuh"
#include "../MHFEM/MHFEM.cuh"
#include "../MHFEM/ReconstructVelocityGrid.cuh"
#include "EdgeToEle.cuh"
#include "Particle.cuh"

namespace cuDFNsys
{
__global__ void ParticleMovementOneTimeStepGPUKernel(unsigned long seed,
                                                     float delta_T_,
                                                     float Dispersion_local,
                                                     cuDFNsys::Particle *P_DEV,
                                                     cuDFNsys::Fracture *Frac_DEV,
                                                     cuDFNsys::EdgeToEle *EdgesSharedEle_DEV,
                                                     cuDFNsys::EleCoor *Coordinate2D_Vec_dev_ptr,
                                                     uint *EleToFracID_ptr,
                                                     float *velocity_ptr,
                                                     int count)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= count)
        return;

    // velocity vector of the grid
    uint EleID = P_DEV[i].ElementID; // from 1
    //uint FracID = EleToFracID_ptr[EleID - 1];                               // from 0
    uint3 EdgeNO = make_uint3(EleID * 3 - 3, EleID * 3 - 2, EleID * 3 - 1); // from 0
    float2 Position_p = P_DEV[i].Position2D;
    float3 Veloc_triangle = make_float3(velocity_ptr[EdgeNO.x],
                                        velocity_ptr[EdgeNO.y],
                                        velocity_ptr[EdgeNO.z]);
    float2 Vertex_Triangle[3];
    Vertex_Triangle[0] = make_float2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[0], Coordinate2D_Vec_dev_ptr[EleID - 1].y[0]);
    Vertex_Triangle[1] = make_float2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[1], Coordinate2D_Vec_dev_ptr[EleID - 1].y[1]);
    Vertex_Triangle[2] = make_float2(Coordinate2D_Vec_dev_ptr[EleID - 1].x[2], Coordinate2D_Vec_dev_ptr[EleID - 1].y[2]);

    float2 Veloc_p = cuDFNsys::ReconstructVelocityGrid(Position_p, Vertex_Triangle, Veloc_triangle);

    // move the particle
    curandState state;

    curand_init(seed, i, 0, &state);

    Position_p.x = Position_p.x + Veloc_p.x * delta_T_ + curand_normal(&state) * sqrt(2.0f * Dispersion_local * delta_T_);
    Position_p.y = Position_p.y + Veloc_p.y * delta_T_ + curand_normal(&state) * sqrt(2.0f * Dispersion_local * delta_T_);

    //printf("%f %f\n", z_1, z_2);

    // check if the new position is in the grid
    bool InGrid = cuDFNsys::IfPntInside2DConvexPoly(Position_p, Vertex_Triangle, 3);
    int edgeNO_local = 0;
    bool OnGridBound = cuDFNsys::IfPntLiesOnBound2DConvexPolyReturnEdgeNO(Position_p, Vertex_Triangle, 3, _TOL_ParticleOnGridBound, &edgeNO_local);

    if (InGrid == true && OnGridBound == false)
    {
        P_DEV[i].Position2D = Position_p;
        // we do not need to change grid ID in this case
        return;
    }

    // check which edge (of the grid) that the connection of last and current positions intersects with?
    // if XX crosses the edge, I will stop the particle on the edge
    // determine the next grid
    int edgeNO_currentPos = 0;
    uint *edgeNO_ptr = &(EdgeNO.x);

    if (InGrid == false && OnGridBound == true) // the new position lies on an edge
    {
        P_DEV[i].Position2D = Position_p;
        edgeNO_currentPos = edgeNO_ptr[edgeNO_local];
    }
    else // the XX crosses an edge
    {
        for (uint ik = 0; ik < 3; ++ik)
        {
            float2 XX[2] = {P_DEV[i].Position2D, Position_p};
            float2 Seg_H[2] = {Vertex_Triangle[ik], Vertex_Triangle[(ik + 1) % 3]};

            float distance1 = cuDFNsys::DistancePnt2DSeg(P_DEV[i].Position2D, Seg_H);
            float distance2 = cuDFNsys::DistancePnt2DSeg(Position_p, Seg_H);

            if (abs(distance1) < _TOL_ParticleOnGridBound && abs(distance2) > _TOL_ParticleOnGridBound)
                continue; // in case of the last position is on an edge of a grid, but the current position is on anther grid
                          // namely, the XX reaches too far

            float2 Intersection[2];
            int sign_o;
            bool Pl = cuDFNsys::IntersectionTwo2DSegs(Seg_H, XX, Intersection, &sign_o, _TOL_ParticleOnGridBound);

            if (Pl == false && ik == 2)
            {
                printf("Warning : I cannot find which edge of a grid that intersects with XX\nThe grid may be a boundary grid, and the movement is out of the fracture\nSo I just cancel this time step\n");
                return;
            }

            if (Pl == true)
            {
                edgeNO_currentPos = edgeNO_ptr[ik];

                if (sign_o == 1)
                    P_DEV[i].Position2D = Intersection[0]; // stop the particle on the edge

                if (sign_o == 2)
                {
                    for (uint j = 0; j < 2; ++j)
                    {
                        float2 A_ = make_float2(P_DEV[i].Position2D.x - Intersection[j].x, P_DEV[i].Position2D.y - Intersection[j].y);
                        float distance_a = sqrt(A_.x * A_.x + A_.y * A_.y);

                        if (distance_a > _TOL_ParticleOnGridBound)
                        {
                            // stop the particle on the edge
                            P_DEV[i].Position2D = Intersection[j];
                        }

                        if (j == 1)
                        {
                            printf("Warning : I cannot find the another end of XX\nSo I just cancel this time step\n");
                            return;
                        }
                    }
                }
                break;
            }
        }
    };
    //----- calculate the elapsed time from last position to current position, which is <= delta_t

    //----- determine where to go next in the remanning time
    uint EleNUMshared = EdgesSharedEle_DEV[edgeNO_currentPos].NumSharedEle;
    if (EleNUMshared == 2)
    {
        P_DEV[i].ElementID = (EdgesSharedEle_DEV[edgeNO_currentPos].EleID[0] == P_DEV[i].ElementID == P_DEV[i].ElementID ? EdgesSharedEle_DEV[edgeNO_currentPos].EleID[1] : EdgesSharedEle_DEV[edgeNO_currentPos].EleID[0]);
        return;
    }

    if (EleNUMshared > 2)
    {
        // flux weights, velocity of grid center
        uint *ElementArray = new uint[EleNUMshared - 1];
        uint jk = 0;

        for (uint j = 0; j < EleNUMshared; ++j)
            if (EdgesSharedEle_DEV[edgeNO_currentPos].EleID[j] != P_DEV[i].ElementID)
            {
                ElementArray[jk] = EdgesSharedEle_DEV[edgeNO_currentPos].EleID[j];
                jk++;
            };

        float *Velocity_grids = new float[EleNUMshared - 1];
        float total_flux = 0;

        for (uint j = 0; j < EleNUMshared - 1; ++j)
        {
            float2 Center_ =
                make_float2(0.5 * (Coordinate2D_Vec_dev_ptr[ElementArray[j] - 1].x[0] + Coordinate2D_Vec_dev_ptr[ElementArray[j] - 1].x[1] + Coordinate2D_Vec_dev_ptr[ElementArray[j] - 1].x[2]),
                            0.5 * (Coordinate2D_Vec_dev_ptr[ElementArray[j] - 1].y[0] + Coordinate2D_Vec_dev_ptr[ElementArray[j] - 1].y[1] + Coordinate2D_Vec_dev_ptr[ElementArray[j] - 1].y[2]));

            float2 Vertex_Triangle_c[3];
            Vertex_Triangle_c[0] = make_float2(Coordinate2D_Vec_dev_ptr[ElementArray[j] - 1].x[0], Coordinate2D_Vec_dev_ptr[ElementArray[j] - 1].y[0]);
            Vertex_Triangle_c[1] = make_float2(Coordinate2D_Vec_dev_ptr[ElementArray[j] - 1].x[1], Coordinate2D_Vec_dev_ptr[ElementArray[j] - 1].y[1]);
            Vertex_Triangle_c[2] = make_float2(Coordinate2D_Vec_dev_ptr[ElementArray[j] - 1].x[2], Coordinate2D_Vec_dev_ptr[ElementArray[j] - 1].y[2]);

            uint3 EdgeNO_c = make_uint3(ElementArray[j] * 3 - 3, ElementArray[j] * 3 - 2, ElementArray[j] * 3 - 1); // from 0

            float3 Veloc_triangle_c = make_float3(velocity_ptr[EdgeNO_c.x],
                                                  velocity_ptr[EdgeNO_c.y],
                                                  velocity_ptr[EdgeNO_c.z]);

            float2 Veloc_ele_c = cuDFNsys::ReconstructVelocityGrid(Center_, Vertex_Triangle_c, Veloc_triangle_c);
            Velocity_grids[j] = sqrt(Veloc_ele_c.x * Veloc_ele_c.x + Veloc_ele_c.y * Veloc_ele_c.y);
            total_flux += Velocity_grids[j];
        };

        float rand_0_1 = curand_uniform(&state);

        float flux_accum_tmp = 0;

        for (uint j = 0; j < EleNUMshared - 1; ++j)
        {
            flux_accum_tmp += Velocity_grids[j];

            float probility_ = flux_accum_tmp / total_flux;

            if (rand_0_1 < probility_)
            {
                P_DEV[i].ElementID = ElementArray[j];
                break;
            }
        };

        delete[] ElementArray;
        ElementArray = NULL;

        delete[] Velocity_grids;
        Velocity_grids = NULL;
        return;
    };

    // this case is the most trick one!
    // the edge is a bound edge! the velocity
    if (EleNUMshared == 1)
    {
        return;
    }
};
}; // namespace cuDFNsys
