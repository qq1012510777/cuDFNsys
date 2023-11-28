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

#include "ParticleTransport/Roate2DPositionTo3D.cuh"

// ====================================================
// NAME:        Roate2DPositionTo3D
// DESCRIPTION: Rotate a 2D particle position to 3D
// AUTHOR:      Tingchang YIN
// DATE:        15/11/2022
// ====================================================
template <typename T>
__host__ __device__ cuDFNsys::Vector3<T>
cuDFNsys::Roate2DPositionTo3D(cuDFNsys::Vector2<T> PositionP,
                              cuDFNsys::Fracture<T> OneFrac)
{
    cuDFNsys::Vector3<T> P_3D =
        cuDFNsys::MakeVector3(PositionP.x, PositionP.y, (T)0.0);

    T RK_2[3][3];
    OneFrac.RoationMatrix(RK_2, 23);

    P_3D = cuDFNsys::ProductSquare3Vector3<T>(RK_2, P_3D);
    P_3D.x += OneFrac.Center.x;
    P_3D.y += OneFrac.Center.y;
    P_3D.z += OneFrac.Center.z;

    return P_3D;
}; // Roate2DPositionTo3D
template __host__ __device__ cuDFNsys::Vector3<double>
cuDFNsys::Roate2DPositionTo3D<double>(cuDFNsys::Vector2<double> PositionP,
                                      cuDFNsys::Fracture<double> OneFrac);
template __host__ __device__ cuDFNsys::Vector3<float>
cuDFNsys::Roate2DPositionTo3D<float>(cuDFNsys::Vector2<float> PositionP,
                                     cuDFNsys::Fracture<float> OneFrac);

// ====================================================
// NAME:        Roate3DPositionTo2D
// DESCRIPTION: Rotate a 3D particle position to 2D
// AUTHOR:      Tingchang YIN
// DATE:        09/11/2023
// ====================================================
template <typename T>
__host__ __device__ cuDFNsys::Vector2<T>
cuDFNsys::Roate3DPositionTo2D(cuDFNsys::Vector3<T> PositionP,
                              cuDFNsys::Fracture<T> OneFrac)
{
    T RK_2[3][3];
    OneFrac.RoationMatrix(RK_2, 32);

    cuDFNsys::Vector3<T> P_3D = PositionP;

    P_3D.x -= OneFrac.Center.x;
    P_3D.y -= OneFrac.Center.y;
    P_3D.z -= OneFrac.Center.z;

    P_3D = cuDFNsys::ProductSquare3Vector3<T>(RK_2, P_3D);

    return cuDFNsys::MakeVector2(P_3D.x, P_3D.y);
};
template __host__ __device__ cuDFNsys::Vector2<double>
cuDFNsys::Roate3DPositionTo2D(cuDFNsys::Vector3<double> PositionP,
                              cuDFNsys::Fracture<double> OneFrac);
template __host__ __device__ cuDFNsys::Vector2<float>
cuDFNsys::Roate3DPositionTo2D(cuDFNsys::Vector3<float> PositionP,
                              cuDFNsys::Fracture<float> OneFrac);

// ====================================================
// NAME:        RotateLineSegToPlane
// DESCRIPTION: RotateLineSegToPlane
// AUTHOR:      Tingchang YIN
// DATE:        09/11/2023
// ====================================================
template <typename T>
__host__ __device__ bool cuDFNsys::RotateLineSegToPlane(
    cuDFNsys::Vector3<T> foot, cuDFNsys::Vector3<T> &end,
    cuDFNsys::Vector3<T> Plane[3], int LocalEdgeOnCEle, T outletcoordinate,
    uint Dir_flow)
{

    T distence_s = cuDFNsys::DistancePnt3DPlane<T>(Plane, foot);
    // printf("distance: %.30f\n", distence_s);
    if (distence_s > 1e-4)
        return false;

    cuDFNsys::Vector3<T> D1_triangle = cuDFNsys::MakeVector3(
        Plane[(LocalEdgeOnCEle + 2) % 3].x - Plane[LocalEdgeOnCEle].x,
        Plane[(LocalEdgeOnCEle + 2) % 3].y - Plane[LocalEdgeOnCEle].y,
        Plane[(LocalEdgeOnCEle + 2) % 3].z - Plane[LocalEdgeOnCEle].z);
    cuDFNsys::Vector3<T> D2_triangle = cuDFNsys::MakeVector3(
        Plane[(LocalEdgeOnCEle + 1) % 3].x - Plane[LocalEdgeOnCEle].x,
        Plane[(LocalEdgeOnCEle + 1) % 3].y - Plane[LocalEdgeOnCEle].y,
        Plane[(LocalEdgeOnCEle + 1) % 3].z - Plane[LocalEdgeOnCEle].z);
    cuDFNsys::Vector3<T> Normal_triangle =
        cuDFNsys::CrossProductVector3<T>(D1_triangle, D2_triangle);
    T norm_1 = pow(Normal_triangle.x * Normal_triangle.x +
                       Normal_triangle.y * Normal_triangle.y +
                       Normal_triangle.z * Normal_triangle.z,
                   0.5);
    Normal_triangle.x /= norm_1;
    Normal_triangle.y /= norm_1;
    Normal_triangle.z /= norm_1;

    cuDFNsys::Vector3<T> D1_trajectory = cuDFNsys::MakeVector3(
        end.x - Plane[LocalEdgeOnCEle].x, end.y - Plane[LocalEdgeOnCEle].y,
        end.z - Plane[LocalEdgeOnCEle].z);
    cuDFNsys::Vector3<T> Normal_trajectory =
        cuDFNsys::CrossProductVector3<T>(D1_trajectory, D2_triangle);
    T norm_2 = pow(Normal_trajectory.x * Normal_trajectory.x +
                       Normal_trajectory.y * Normal_trajectory.y +
                       Normal_trajectory.z * Normal_trajectory.z,
                   0.5);
    Normal_trajectory.x /= norm_2;
    Normal_trajectory.y /= norm_2;
    Normal_trajectory.z /= norm_2;

    T angle_s = acos(Normal_triangle.x * Normal_trajectory.x +
                     Normal_triangle.y * Normal_trajectory.y +
                     Normal_triangle.z * Normal_trajectory.z);

    // printf("angle_s: %.30f, norm of trajectory: %.30f, norm_triangle: %.30f, "
    //        "%.30f %.30f; Normal_trajectory: %.30f, %.30f, %.30f; dot_product: "
    //        "%.30f\n\n",
    //        angle_s,
    //        pow((foot.x - end.x) * (foot.x - end.x) +
    //                (foot.y - end.y) * (foot.y - end.y) +
    //                (foot.z - end.z) * (foot.z - end.z),
    //            0.5),
    //        Normal_triangle.x, Normal_triangle.y, Normal_triangle.z,
    //        Normal_trajectory.x, Normal_trajectory.y, Normal_trajectory.z,
    //        Normal_triangle.x * Normal_trajectory.x +
    //            Normal_triangle.y * Normal_trajectory.y +
    //            Normal_triangle.z * Normal_trajectory.z);

    if (angle_s != angle_s)
    {
        if ((abs(Normal_triangle.x - 1) < 1e-5 &&
             abs(Normal_triangle.y) < 1e-5 && abs(Normal_triangle.z) < 1e-5) ||
            (abs(Normal_triangle.x) < 1e-5 &&
             abs(Normal_triangle.y - 1) < 1e-5 &&
             abs(Normal_triangle.z) < 1e-5) ||
            (abs(Normal_triangle.x) < 1e-5 && abs(Normal_triangle.y) < 1e-5 &&
             abs(Normal_triangle.z - 1) < 1e-5))
            angle_s = 0.;
        else
        {
            printf("The angle is NaN: %.30f, norm_1: %.30f, Normal_triangle: "
                   "%.30f, %.30f, %.30f\n",
                   angle_s, norm_1, Normal_triangle.x, Normal_triangle.y,
                   Normal_triangle.z);
            return false;
        }
    }

    cuDFNsys::Vector3<T> Tar = cuDFNsys::MakeVector3(
        end.x - Plane[LocalEdgeOnCEle].x, end.y - Plane[LocalEdgeOnCEle].y,
        end.z - Plane[LocalEdgeOnCEle].z);
    cuDFNsys::Vector3<T> RotateAxis = cuDFNsys::MakeVector3(
        Plane[(LocalEdgeOnCEle + 1) % 3].x - Plane[LocalEdgeOnCEle].x,
        Plane[(LocalEdgeOnCEle + 1) % 3].y - Plane[LocalEdgeOnCEle].y,
        Plane[(LocalEdgeOnCEle + 1) % 3].z - Plane[LocalEdgeOnCEle].z);

    T norm_rot = pow(RotateAxis.x * RotateAxis.x + RotateAxis.y * RotateAxis.y +
                         RotateAxis.z * RotateAxis.z,
                     0.5);
    RotateAxis.x /= norm_rot;
    RotateAxis.y /= norm_rot;
    RotateAxis.z /= norm_rot;

    //printf("angle_s: %.30f\n", angle_s);
    bool If_found = false;
    for (int i_s = 0; i_s < 2; ++i_s)
    {
        T angle_fsf = angle_s;
        if (i_s == 1)
            angle_fsf = M_PI - angle_fsf;

        for (int i = 0; i < 5; ++i)
        {
            cuDFNsys::Quaternion<T> qua;
            T angle_IO = angle_fsf + M_PI * i;

            qua = qua.DescribeRotation(RotateAxis, angle_IO);
            cuDFNsys::Vector3<T> tar1 = qua.Rotate(Tar);

            tar1.x += Plane[LocalEdgeOnCEle].x;
            tar1.y += Plane[LocalEdgeOnCEle].y;
            tar1.z += Plane[LocalEdgeOnCEle].z;

            T distence_f = cuDFNsys::DistancePnt3DPlane<T>(Plane, tar1);
            // printf("i_s: %d, i: %d, distence_f: %.30f, tar1: %.30f, %.30f, %.30f\nTar: %.30f, "
            //        "%.30f, %.30f\nPlane[LocalEdgeOnCEle]: %.30f, "
            //        "%.30f, %.30f\n\n", i_s, i,
            //      distence_f, tar1.x, tar1.y, tar1.z,
            //      Tar.x + Plane[LocalEdgeOnCEle].x,
            //      Tar.y + Plane[LocalEdgeOnCEle].y,
            //      Tar.z + Plane[LocalEdgeOnCEle].z, Plane[LocalEdgeOnCEle].x,
            //      Plane[LocalEdgeOnCEle].y, Plane[LocalEdgeOnCEle].z);
            if (distence_f < 1e-5)
            {
                T *tmp_t = &(tar1.x);
                if (abs(tmp_t[Dir_flow]) < abs(outletcoordinate))
                {
                    Tar = tar1;
                    If_found = true;
                    break;
                }
            }
            if (i == 4 && i_s == 1)
                return false;
        }
        if (If_found)
            break;
    }

    //printf("foot: %.30f, %.30f, %.30f\n", foot.x, foot.y, foot.z);
    //printf("end: %.30f, %.30f, %.30f\n", end.x, end.y, end.z);
    //printf("Triangle 1: %.30f, %.30f, %.30f\n", Plane[0].x, Plane[0].y,
    //       Plane[0].z);
    //printf("Triangle 2: %.30f, %.30f, %.30f\n", Plane[1].x, Plane[1].y,
    //       Plane[1].z);
    //printf("Triangle 3: %.30f, %.30f, %.30f\n", Plane[2].x, Plane[2].y,
    //       Plane[2].z);
    //printf("RotateAxis 1: %.30f, %.30f, %.30f\n", Plane[LocalEdgeOnCEle].x,
    //       Plane[LocalEdgeOnCEle].y, Plane[LocalEdgeOnCEle].z);
    //printf("RotateAxis 2: %.30f, %.30f, %.30f\n",
    //       Plane[(LocalEdgeOnCEle + 1) % 3].x,
    //       Plane[(LocalEdgeOnCEle + 1) % 3].y,
    //       Plane[(LocalEdgeOnCEle + 1) % 3].z);
    //printf("norm 1: %.30f, %.30f, %.30f\n", Normal_triangle.x,
    //       Normal_triangle.y, Normal_triangle.z);
    //printf("norm 2: %.30f, %.30f, %.30f\n", Normal_trajectory.x,
    //       Normal_trajectory.y, Normal_trajectory.z);
    //printf("angle_degree: %.30f\n", angle_s * 180.0 / M_PI);
    //printf("Tar: %.30f, %.30f, %.30f\n", Tar.x, Tar.y, Tar.z);
    end = Tar;
    return true;
};
template __host__ __device__ bool cuDFNsys::RotateLineSegToPlane<double>(
    cuDFNsys::Vector3<double> foot, cuDFNsys::Vector3<double> &end,
    cuDFNsys::Vector3<double> Plane[3], int LocalEdgeOnCEle,
    double outletcoordinate, uint Dir_flow);
template __host__ __device__ bool cuDFNsys::RotateLineSegToPlane<float>(
    cuDFNsys::Vector3<float> foot, cuDFNsys::Vector3<float> &end,
    cuDFNsys::Vector3<float> Plane[3], int LocalEdgeOnCEle,
    float outletcoordinate, uint Dir_flow);