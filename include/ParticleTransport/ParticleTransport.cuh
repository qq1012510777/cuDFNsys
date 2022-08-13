///////////////////////////////////////////////////////////////////
// NAME:              ParticleTransport.cuh
//
// PURPOSE:           perform particle tracking with random walk method
//
// FUNCTIONS/OBJECTS: ParticleTransport
//
// AUTHOR:            Tingchang YIN
///////////////////////////////////////////////////////////////////
#pragma once
#include "EdgeToEle.cuh"
#include "NeighborEle.cuh"
#include "Particle.cuh"
#include "ParticleMovementOneTimeStepGPUKernel.cuh"
#include "PredicateNumOfReachedOutletParticles.cuh"
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

namespace cuDFNsys
{
template <typename T>
class ParticleTransport
{
public:
    int NumParticles = 0;
    thrust::host_vector<cuDFNsys::Particle<T>> ParticlePlumes;
    thrust::host_vector<cuDFNsys::EdgeToEle> EdgesSharedEle;
    thrust::host_vector<cuDFNsys::NeighborEle> NeighborEleOfOneEle;

    uint Dir = 2;

private:
    string ParticlePosition = "ParticlePositionResult/ParticlePosition";

public:
    ParticleTransport(unsigned long seed,
                      const int &NumOfParticles,
                      const int &NumTimeStep,
                      T delta_T_,
                      T Dispersion_local,
                      thrust::host_vector<cuDFNsys::Fracture<T>> Fracs,
                      cuDFNsys::Mesh<T> mesh,
                      const cuDFNsys::MHFEM<T> &fem,
                      uint Dir_flow,
                      T outletcoordinate,
                      const string &Particle_mode,
                      const string &Injection_mode)
    {
        this->Dir = Dir_flow;

        this->IdentifyEdgesSharedEle(mesh);

        this->IdentifyNeighborElements(mesh);

        this->InitilizeParticles(NumOfParticles, mesh, fem, Injection_mode);

        this->OutputParticleInfoStepByStep(ParticlePosition, 0,
                                           delta_T_,
                                           Dispersion_local,
                                           Particle_mode,
                                           Injection_mode,
                                           Fracs, mesh);
        //cout << NumTimeStep << ", " << delta_T_ << ", " << Dispersion_local << " ______ \n";
        this->ParticleMovement(seed, 1, NumTimeStep, delta_T_, Dispersion_local, Particle_mode,
                               Injection_mode, Fracs, mesh, fem, outletcoordinate);
    };

    ParticleTransport(unsigned long seed,
                      const int &NumTimeStep,
                      thrust::host_vector<cuDFNsys::Fracture<T>> Fracs,
                      cuDFNsys::Mesh<T> mesh,
                      const cuDFNsys::MHFEM<T> &fem,
                      uint Dir_flow,
                      T outletcoordinate)
    {
        cuDFNsys::MatlabAPI M1;

        this->Dir = Dir_flow;

        this->IdentifyEdgesSharedEle(mesh);

        this->IdentifyNeighborElements(mesh);

        MatrixXd PreNumSteps;
        string matfile = ParticlePosition + "_step_" + cuDFNsys::ToStringWithWidth(0, 7) + ".mat";
        M1.ReadMat(matfile, "NumOfSteps", PreNumSteps);
        uint ExistingNumsteps = (uint)PreNumSteps(0, 0);

        string Particle_mode = M1.ReadMatString(matfile, "Particle_mode");
        string Injection_mode = M1.ReadMatString(matfile, "Injection_mode");

        MatrixXd LastStepParticle;
        M1.ReadMat(ParticlePosition + "_step_" + cuDFNsys::ToStringWithWidth(ExistingNumsteps, 7) + ".mat", "particle_position_3D_step_" + std::to_string(ExistingNumsteps),
                   LastStepParticle);

        this->NumParticles = LastStepParticle.rows();
        this->ParticlePlumes.resize(this->NumParticles);
        for (uint i = 0; i < this->NumParticles; ++i)
        {
            this->ParticlePlumes[i].Position2D.x = LastStepParticle(i, 5);
            this->ParticlePlumes[i].Position2D.y = LastStepParticle(i, 6);
            this->ParticlePlumes[i].ElementID = LastStepParticle(i, 4);
            this->ParticlePlumes[i].IfReachOutletPlane = (LastStepParticle(i, 3) == 0 ? false : true);
        };

        MatrixXd delta_T_b, Dispersion_local_b;
        M1.ReadMat(matfile, "Delta_T",
                   delta_T_b);
        M1.ReadMat(matfile, "Dispersion_local",
                   Dispersion_local_b);

        T delta_T_ = (T)delta_T_b(0, 0),
          Dispersion_local = (T)Dispersion_local_b(0, 0);
        
        this->ParticleMovement(seed, ExistingNumsteps + 1, NumTimeStep,
                               delta_T_, Dispersion_local,
                               Particle_mode,
                               Injection_mode,
                               Fracs, mesh, fem, outletcoordinate);
    };

    void ParticleMovement(unsigned long seed,
                          const int &init_NO_STEP,
                          const int &NumTimeStep,
                          T delta_T_,
                          T Dispersion_local,
                          const string &Particle_mode,
                          const string &Injection_mode,
                          thrust::host_vector<cuDFNsys::Fracture<T>> Fracs,
                          cuDFNsys::Mesh<T> mesh,
                          const cuDFNsys::MHFEM<T> &fem,
                          T outletcoordinate)
    {

        thrust::device_vector<cuDFNsys::Particle<T>> ParticlePlumes_DEV = this->ParticlePlumes;
        thrust::device_vector<cuDFNsys::Fracture<T>> Fracsvec_DEV = Fracs;
        thrust::device_vector<cuDFNsys::EdgeToEle> EdgesSharedEle_vec_dev = this->EdgesSharedEle;
        thrust::device_vector<uint> ElementFracTag_DEV = mesh.ElementFracTag;
        thrust::device_vector<T> Velocity_sep_edge;
        Velocity_sep_edge.reserve(fem.VelocityNormalScalarSepEdges.rows());
        thrust::device_vector<cuDFNsys::EleCoor<T>> Coordinate2D_Vec_dev = mesh.Coordinate2D;
        thrust::device_vector<cuDFNsys::NeighborEle> NeighborEleOfOneEle_dev = this->NeighborEleOfOneEle;

        Eigen::MatrixXd Vg = fem.VelocityNormalScalarSepEdges;
        double *vc = Vg.data();
        Velocity_sep_edge.insert(Velocity_sep_edge.end(), &vc[0], &vc[fem.VelocityNormalScalarSepEdges.rows()]);

        cuDFNsys::Particle<T> *P_DEV = thrust::raw_pointer_cast(ParticlePlumes_DEV.data());
        cuDFNsys::Fracture<T> *Frac_DEV = thrust::raw_pointer_cast(Fracsvec_DEV.data());
        cuDFNsys::EdgeToEle *EdgesSharedEle_DEV = thrust::raw_pointer_cast(EdgesSharedEle_vec_dev.data());
        uint *EleToFracID_ptr = thrust::raw_pointer_cast(ElementFracTag_DEV.data());
        T *velocity_ptr = thrust::raw_pointer_cast(Velocity_sep_edge.data());
        cuDFNsys::EleCoor<T> *Coordinate2D_Vec_dev_ptr = thrust::raw_pointer_cast(Coordinate2D_Vec_dev.data());
        cuDFNsys::NeighborEle *NeighborEleOfOneEle_dev_ptr = thrust::raw_pointer_cast(NeighborEleOfOneEle_dev.data());

        for (uint i = init_NO_STEP; i <= NumTimeStep + init_NO_STEP; ++i)
        {
            thrust::host_vector<uint> Particle_runtime_error(this->ParticlePlumes.size(), 0);
            thrust::device_vector<uint> Particle_runtime_error_dev = Particle_runtime_error;
            uint *Particle_runtime_error_dev_pnt = thrust::raw_pointer_cast(Particle_runtime_error_dev.data());

            time_t t;
            time(&t);

            cout << "The Step " << i << endl;
            double istart_b = cuDFNsys::CPUSecond();
            ParticleMovementOneTimeStepGPUKernel<T><<<this->NumParticles / 256 + 1, 256>>>((unsigned long)t + (unsigned long)(i * 100),
                                                                                           delta_T_,
                                                                                           Dispersion_local,
                                                                                           P_DEV,
                                                                                           Frac_DEV,
                                                                                           EdgesSharedEle_DEV,
                                                                                           Coordinate2D_Vec_dev_ptr,
                                                                                           NeighborEleOfOneEle_dev_ptr,
                                                                                           EleToFracID_ptr,
                                                                                           velocity_ptr,
                                                                                           this->Dir,
                                                                                           outletcoordinate,
                                                                                           NumParticles,
                                                                                           mesh.Element3D.size(),
                                                                                           i,
                                                                                           Particle_runtime_error_dev_pnt);
            cudaDeviceSynchronize();
            Particle_runtime_error = Particle_runtime_error_dev;

            double ielaps_b = cuDFNsys::CPUSecond() - istart_b;

            uint sum = thrust::reduce(Particle_runtime_error.begin(), Particle_runtime_error.end(), (uint)0, thrust::plus<uint>());
            if (sum != 0)
                throw cuDFNsys::ExceptionsPause("Error happens when particle is moving!\n");

            this->ParticlePlumes = ParticlePlumes_DEV;
            double istart = cuDFNsys::CPUSecond();

            int NumReachedParticle = thrust::count_if(this->ParticlePlumes.begin(),
                                                      this->ParticlePlumes.end(),
                                                      cuDFNsys::PredicateNumOfReachedOutletParticles<T>());
            double ielaps = cuDFNsys::CPUSecond() - istart;
            cout << NumReachedParticle << "/" << this->ParticlePlumes.size() << " reached outlet plane, running time: " << ielaps_b << "; counting time: " << ielaps << "s\n";
            this->OutputParticleInfoStepByStep(ParticlePosition, i,
                                               delta_T_,
                                               Dispersion_local,
                                               Particle_mode,
                                               Injection_mode,
                                               Fracs, mesh);

            cout << "\n";
        }
    };

    void OutputParticleInfoStepByStep(const string &mat_key_r,
                                      const uint &StepNO,
                                      const T delta_T,
                                      const T Dispersion_local,
                                      const string &Particle_mode,
                                      const string &Injection_mode,
                                      thrust::host_vector<cuDFNsys::Fracture<T>> Fracs,
                                      cuDFNsys::Mesh<T> mesh)
    {
        cuDFNsys::MatlabAPI M1;

        string mat_key = mat_key_r + "_step_" + cuDFNsys::ToStringWithWidth(StepNO, 7) + ".mat";

        if (StepNO == 0)
        {
            M1.WriteMatString(mat_key, "w", Particle_mode, "Particle_mode");
            M1.WriteMatString(mat_key, "u", Injection_mode, "Injection_mode");
            T po[1] = {delta_T};
            M1.WriteMat(mat_key, "u", 1,
                        1, 1, po, "Delta_T");
            po[0] = {Dispersion_local};
            M1.WriteMat(mat_key, "u", 1,
                        1, 1, po, "Dispersion_local");
        }

        uint Step[1] = {StepNO};
        M1.WriteMat(mat_key_r + "_step_" + cuDFNsys::ToStringWithWidth(0, 7) + ".mat", "u", 1,
                    1, 1, Step, "NumOfSteps");

        T *particle_position_3D;
        particle_position_3D = new T[this->NumParticles * 7];
        if (particle_position_3D == NULL)
        {
            string AS = "Alloc error in ParticleTransport::OutputParticleInfoStepByStep\n";
            throw cuDFNsys::ExceptionsPause(AS);
        }

        for (size_t i = 0; i < this->NumParticles; ++i)
        {
            cuDFNsys::Vector3<T> tmpPos = cuDFNsys::MakeVector3(this->ParticlePlumes[i].Position2D.x,
                                                                this->ParticlePlumes[i].Position2D.y,
                                                                (T)0);

            uint elementID = this->ParticlePlumes[i].ElementID; // from 1
            uint FracID = mesh.ElementFracTag[elementID - 1];   // from 0

            T Rotate2DTo3D[3][3];
            Fracs[FracID].RoationMatrix(Rotate2DTo3D, 23);

            tmpPos = cuDFNsys::ProductSquare3Vector3<T>(Rotate2DTo3D, tmpPos);

            tmpPos = cuDFNsys::MakeVector3(tmpPos.x + Fracs[FracID].Center.x,
                                           tmpPos.y + Fracs[FracID].Center.y,
                                           tmpPos.z + Fracs[FracID].Center.z);

            particle_position_3D[i] = tmpPos.x;
            particle_position_3D[i + this->NumParticles] = tmpPos.y;
            particle_position_3D[i + this->NumParticles * 2] = tmpPos.z;
            particle_position_3D[i + this->NumParticles * 3] = (this->ParticlePlumes[i].IfReachOutletPlane == false ? 0 : 1);
            particle_position_3D[i + this->NumParticles * 4] = this->ParticlePlumes[i].ElementID;
            particle_position_3D[i + this->NumParticles * 5] = this->ParticlePlumes[i].Position2D.x;
            particle_position_3D[i + this->NumParticles * 6] = this->ParticlePlumes[i].Position2D.y;
        }

        string mode = "w";

        if (StepNO == 0)
            mode = "u";

        M1.WriteMat(mat_key, mode, this->NumParticles * 7,
                    this->NumParticles, 7, particle_position_3D,
                    "particle_position_3D_step_" + to_string(StepNO));
        delete[] particle_position_3D;
        particle_position_3D = NULL;
    };

    void MatlabPlot(const string &mat_key,
                    const string &command_key,
                    const cuDFNsys::Mesh<T> &mesh,
                    const cuDFNsys::MHFEM<T> &fem,
                    const T &L)
    {
        // cuDFNsys::MatlabAPI M1;
        //
        // size_t node_num = mesh.Coordinate3D.size();
        //
        // T *ptr_coordinates_3D;
        // ptr_coordinates_3D = new T[node_num * 3];
        //
        // if (ptr_coordinates_3D == NULL)
        // {
        //     string AS = "Alloc error in ParticleTransport::MatlabPlot\n";
        //     throw cuDFNsys::ExceptionsPause(AS);
        // }
        //
        // for (size_t i = 0; i < node_num; ++i)
        // {
        //     ptr_coordinates_3D[i] = mesh.Coordinate3D[i].x;
        //     ptr_coordinates_3D[i + node_num] = mesh.Coordinate3D[i].y;
        //     ptr_coordinates_3D[i + node_num * 2] = mesh.Coordinate3D[i].z;
        // }
        // M1.WriteMat(mat_key, "w", node_num * 3,
        //             node_num, 3, ptr_coordinates_3D, "coordinate_3D");
        //
        // delete[] ptr_coordinates_3D;
        // ptr_coordinates_3D = NULL;
        //
        // //---------------------
        // size_t ele_num = mesh.Element3D.size();
        // T *ptr_element_3D = new T[ele_num * 3];
        // if (ptr_element_3D == NULL)
        // {
        //     string AS = "Alloc error in ParticleTransport::MatlabPlot\n";
        //     throw cuDFNsys::ExceptionsPause(AS);
        // }
        //
        // for (size_t i = 0; i < ele_num; ++i)
        // {
        //     ptr_element_3D[i] = mesh.Element3D[i].x;
        //     ptr_element_3D[i + ele_num] = mesh.Element3D[i].y;
        //     ptr_element_3D[i + ele_num * 2] = mesh.Element3D[i].z;
        // }
        // M1.WriteMat(mat_key, "u", ele_num * 3,
        //             ele_num, 3, ptr_element_3D, "element_3D");
        //
        // delete[] ptr_element_3D;
        // ptr_element_3D = NULL;
        // //--------------------------
        //
        // T *pressure_ELEs = new T[ele_num];
        // if (pressure_ELEs == NULL)
        // {
        //     string AS = "Alloc error in ParticleTransport::MatlabPlot\n";
        //     throw cuDFNsys::ExceptionsPause(AS);
        // }
        //
        // std::copy(fem.PressureEles.data(),
        //           fem.PressureEles.data() + ele_num,
        //           pressure_ELEs);
        //
        // M1.WriteMat(mat_key, "u", ele_num,
        //             ele_num, 1, pressure_ELEs,
        //             "pressure_eles");
        //
        // delete[] pressure_ELEs;
        // pressure_ELEs = NULL;
        //----------------------------------

        std::ofstream oss(command_key, ios::out);
        oss << "clc;\nclose all;\nclear all;\ncurrentPath = fileparts(mfilename('fullpath'));\n";
        oss << "load('" << mat_key << "');\n";
        oss << "L = 0.5 * " << L << ";\n";
        oss << "cube_frame = [-L, -L, L; -L, L, L; L, L, L; L -L, L; -L, -L, -L; -L, L, -L; L, L, -L; L -L, -L; -L, L, L; -L, L, -L; -L, -L, -L; -L, -L, L; L, L, L; L, L, -L; L, -L, -L; L, -L, L; L, -L, L; L, -L, -L; -L, -L, -L; -L, -L, L; L, L, L; L, L,-L; -L, L, -L; -L,L, L];\n";
        oss << "figure(1); view(3); title('DFN flow (mhfem)'); xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); hold on\n";
        oss << "patch('Vertices', cube_frame, 'Faces', [1, 2, 3, 4;5 6 7 8;9 10 11 12; 13 14 15 16], 'FaceVertexCData', zeros(size(cube_frame, 1), 1), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0); hold on\n";
        oss << endl;
        oss << "patch('Vertices', coordinate_3D, 'Faces', element_3D, 'FaceVertexCData', pressure_eles, 'FaceColor', 'flat', 'EdgeAlpha', 0.2, 'facealpha', 0.9); colorbar; view(3); hold on\n";
        oss << "caxis([" << fem.OutletP << ", " << fem.InletP << "]);\n";
        oss << "xlim([-(0.1 * L + L), (0.1 * L + L)])\nylim([ -(0.1 * L + L), (0.1 * L + L) ])\nzlim([ -(0.1 * L + L), (0.1 * L + L) ]);hold on\n\n";
        oss << "S = load([currentPath, '/" << ParticlePosition << "_step_" << cuDFNsys::ToStringWithWidth(0, 7) << ".mat']);\n";
        oss << "N_steps = S.NumOfSteps;\n";
        oss << "N_particles = size(S.particle_position_3D_step_0, 1);\n";
        oss << "Matrx3D_pso = zeros(N_particles, size(S.particle_position_3D_step_0, 2), N_steps + 1);\nclear S;\n";
        oss << "for i = 0:N_steps\n";
        //oss << "\ttitle(['DFN flow (mhfem); step NO = ', num2str(i)]);\n";
        oss << "\tS = load([currentPath, '/" << ParticlePosition << "_step_', num2str(i, '%07d'), '.mat']);\n";
        oss << "\teval(['Matrx3D_pso(:, :, i + 1) = S.particle_position_3D_step_', num2str(i), '(:, :);']);\n";
        //oss << "\teval(['p_s = scatter3(S.particle_position_3D_step_', num2str(i), '(:, 1), S.particle_position_3D_step_', num2str(i), '(:, 2), S.particle_position_3D_step_', num2str(i), '(:, 3), ''k'', ''o'', ''filled'');']);\n";
        //oss << "\tif (i == 0); pause; else; pause(0.01); end;\n";
        oss << "\tclear S;\n";
        //oss << "\tif (i ~= S.NumOfSteps); delete(p_s); end\n";
        oss << "end\n";

        oss << "ReachedParticleNO = find(Matrx3D_pso(:, 4, N_steps) == 0);\n";
        oss << "AK_1 = [reshape(Matrx3D_pso(:, 1, :), size(Matrx3D_pso, 1), size(Matrx3D_pso, 3))]';\n";
        oss << "AK_2 = [reshape(Matrx3D_pso(:, 2, :), size(Matrx3D_pso, 1), size(Matrx3D_pso, 3))]';\n";
        oss << "AK_3 = [reshape(Matrx3D_pso(:, 3, :), size(Matrx3D_pso, 1), size(Matrx3D_pso, 3))]';\n";
        oss << "AK_1(:, ReachedParticleNO) = [];\n";
        oss << "AK_2(:, ReachedParticleNO) = [];\n";
        oss << "AK_3(:, ReachedParticleNO) = [];\n";
        oss << "newcolors = rand(size(ReachedParticleNO, 1), 3);\n";
        oss << "colororder(newcolors)\n";
        oss << "plot3(AK_1, AK_2, AK_3, 'linewidth', 2);\n\n";

        oss << "figure(2); view(3); title('DFN flow (mhfem)'); xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); hold on\n";
        oss << "patch('Vertices', cube_frame, 'Faces', [1, 2, 3, 4;5 6 7 8;9 10 11 12; 13 14 15 16], 'FaceVertexCData', zeros(size(cube_frame, 1), 1), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0); hold on\n";
        oss << endl;
        oss << "patch('Vertices', coordinate_3D, 'Faces', element_3D, 'FaceVertexCData', pressure_eles, 'FaceColor', 'flat', 'EdgeAlpha', 0.2, 'facealpha', 0.9); colorbar; view(3); hold on\n";
        oss << "caxis([" << fem.OutletP << ", " << fem.InletP << "]);\n";
        oss << "xlim([-(0.1 * L + L), (0.1 * L + L)])\nylim([ -(0.1 * L + L), (0.1 * L + L) ])\nzlim([ -(0.1 * L + L), (0.1 * L + L) ]);hold on\n\n";
        oss << "figure(2)\nfor i = 0:N_steps\n";
        oss << "\ttitle(['DFN flow (mhfem); step NO = ', num2str(i)]);\n";
        oss << "\tp_s = scatter3(Matrx3D_pso(:, 1, i + 1), Matrx3D_pso(:, 2, i + 1), Matrx3D_pso(:, 3, i + 1), 'k', 'o', 'filled');\n";
        oss << "\tif (i == 0); pause; else; pause(0.01); end;\n";
        oss << "\tif (i ~= N_steps); delete(p_s); end\n";
        oss << "end\n";
    };

private:
    void IdentifyEdgesSharedEle(cuDFNsys::Mesh<T> mesh)
    {
        this->EdgesSharedEle.resize(mesh.Element3D.size() * 3);

        // create unordered_map for each edges
        typedef struct SharedElements
        {
            uint ElementID[_NumOfSharedEleAtMost] = {0};
            uint LocalEdgeNO[_NumOfSharedEleAtMost] = {0};
        } SharedEle;

        std::unordered_map<pair<size_t, size_t>, SharedEle, PairHash> UMapEdge2Ele;
        UMapEdge2Ele.reserve(mesh.Element3D.size() * 3);

        for (size_t i = 0; i < mesh.Element3D.size(); ++i)
        {
            uint *pts = &mesh.Element3D[i].x;

            for (size_t j = 0; j < 3; ++j)
            {
                size_t node1 = pts[j];
                size_t node2 = pts[(j + 1) % 3];

                pair<size_t, size_t> key_;
                key_.first = node1 < node2 ? node1 : node2;
                key_.second = node1 > node2 ? node1 : node2;

                if (UMapEdge2Ele.find(key_) == UMapEdge2Ele.end())
                {
                    UMapEdge2Ele[key_].ElementID[0] = i + 1;
                    UMapEdge2Ele[key_].LocalEdgeNO[0] = j;
                }
                else
                {
                    for (size_t k = 0; k < _NumOfSharedEleAtMost; ++k)
                    {
                        if (UMapEdge2Ele[key_].ElementID[k] == 0)
                        {
                            UMapEdge2Ele[key_].ElementID[k] = i + 1;
                            UMapEdge2Ele[key_].LocalEdgeNO[k] = j;
                            break;
                        };

                        if (k == _NumOfSharedEleAtMost - 1)
                        {
                            throw cuDFNsys::ExceptionsIgnore("Find intersections shared by more than four elements!");
                        }
                    }
                };
            }
        };

        // identify
        for (size_t i = 0; i < this->EdgesSharedEle.size(); ++i)
        {
            uint *pts = &(mesh.Element3D[i / 3].x);

            size_t node1 = pts[i % 3];
            size_t node2 = pts[(i % 3 + 1) % 3];

            pair<size_t, size_t> key_;
            key_.first = node1 < node2 ? node1 : node2;
            key_.second = node1 > node2 ? node1 : node2;

            for (size_t k = 0; k < _NumOfSharedEleAtMost; ++k)
            {
                if (UMapEdge2Ele[key_].ElementID[k] != 0)
                {
                    this->EdgesSharedEle[i].EleID[k] = UMapEdge2Ele[key_].ElementID[k];
                    this->EdgesSharedEle[i].LocalEdgeNO[k] = UMapEdge2Ele[key_].LocalEdgeNO[k];
                    this->EdgesSharedEle[i].NumSharedEle++;
                }
            }

            // for (size_t k = 0; k < this->EdgesSharedEle[i].NumSharedEle; ++k)
            // {
            //     cout << UMapEdge2Ele[key_].ElementID[k] << (k == this->EdgesSharedEle[i].NumSharedEle - 1 ? "\n" : ", ");
            //     if(UMapEdge2Ele[key_].ElementID[k] > mesh.Element3D.size())
            //     {
            //         throw cuDFNsys::ExceptionsIgnore("Out of range of element ID\n");
            //     }
            // }
        }
    };

    void InitilizeParticles(const int &NumOfParticles,
                            cuDFNsys::Mesh<T> mesh,
                            const cuDFNsys::MHFEM<T> &fem,
                            const string &Injection_mode)
    {
        thrust::host_vector<uint> NumParticlesEachEdge(mesh.InletEdgeNOLen.size(), 0);

        if (Injection_mode == "Flux-weighted" || Injection_mode == "Resident")
            for (size_t i = 0; i < mesh.InletEdgeNOLen.size(); ++i)
            {
                uint EdgeNO = (uint)mesh.InletEdgeNOLen[i].x;
                T length_ = mesh.InletEdgeNOLen[i].y;

                if (Injection_mode == "Flux-weighted")
                {

                    T proportion_ = abs(fem.VelocityNormalScalarSepEdges(EdgeNO - 1, 0)) * length_ / fem.QIn;
                    /// flux-weighted

                    NumParticlesEachEdge[i] = (uint)round(proportion_ * NumOfParticles);
                }
                else if (Injection_mode == "Resident")
                {
                    // uniform
                    T proportion_ = abs(length_ / fem.InletLength);
                    NumParticlesEachEdge[i] = (uint)round(proportion_ * NumOfParticles);
                }
            }
        else if (Injection_mode == "Point")
        {
            int a;
            srand((unsigned)time(NULL));
            //cout << "NumOfParticles: " << NumOfParticles << endl;
            for (uint i = 0; i < 100; ++i)
            {
                a = rand() % (int)mesh.InletEdgeNOLen.size();
                //printf("a = %d\n", a);
                uint EdgeNO = (uint)mesh.InletEdgeNOLen[a].x;
                if (fem.VelocityNormalScalarSepEdges(EdgeNO - 1, 0) < 0)
                    NumParticlesEachEdge[a] = NumOfParticles;
                break;
            }
        }
        else
            throw cuDFNsys::ExceptionsPause("Undefined particle injection mode!\n");

        this->NumParticles = 0;
        for (size_t i = 0; i < NumParticlesEachEdge.size(); ++i)
        {
            this->NumParticles += NumParticlesEachEdge[i];
        }
        this->ParticlePlumes.resize(NumParticles);

        uint TmpCountParticle = 0;
        for (size_t i = 0; i < NumParticlesEachEdge.size(); ++i)
        {
            uint EdgeNO = (uint)mesh.InletEdgeNOLen[i].x;     // from 1
            uint elementID = (EdgeNO - 1) / 3 + 1;            // from 1
            uint FracID = mesh.ElementFracTag[elementID - 1]; // from 0

            uint LocalEdgeNO = (EdgeNO - 1) % 3;

            // cuDFNsys::Vector2<T> position_ = cuDFNsys::MakeVector2((T)0.5f * (mesh.Coordinate2D[elementID - 1].x[LocalEdgeNO] + mesh.Coordinate2D[elementID - 1].x[(LocalEdgeNO + 1) % 3]),
            //                                                        (T)0.5f * (mesh.Coordinate2D[elementID - 1].y[LocalEdgeNO] + mesh.Coordinate2D[elementID - 1].y[(LocalEdgeNO + 1) % 3]));
            // // move the position a little toward to the interior of grid
            // {
            //     cuDFNsys::Vector2<T> inwardNormal = cuDFNsys::MakeVector2(-(mesh.Coordinate2D[elementID - 1].y[(LocalEdgeNO + 1) % 3] - mesh.Coordinate2D[elementID - 1].y[LocalEdgeNO]),
            //                                                               mesh.Coordinate2D[elementID - 1].x[(LocalEdgeNO + 1) % 3] - mesh.Coordinate2D[elementID - 1].x[LocalEdgeNO]);
            //     T norm_inward = sqrt(inwardNormal.x * inwardNormal.x + inwardNormal.y * inwardNormal.y);
            //     inwardNormal.x /= norm_inward;
            //     inwardNormal.y /= norm_inward;
            //     position_.x += (inwardNormal.x * 1e-5);
            //     position_.y += (inwardNormal.y * 1e-5);
            // };
            // cuDFNsys::Particle<T> tmpP;
            // tmpP.ElementID = elementID;
            // tmpP.Position2D = position_;
            // thrust::host_vector<cuDFNsys::Particle<T>> tmpPVEC(NumParticlesEachEdge[i], tmpP);

            thrust::host_vector<cuDFNsys::Particle<T>> tmpPVEC(NumParticlesEachEdge[i]);
            cuDFNsys::Vector2<T> Vec2D;
            Vec2D.x = mesh.Coordinate2D[elementID - 1].x[(LocalEdgeNO + 1) % 3] - mesh.Coordinate2D[elementID - 1].x[LocalEdgeNO];
            Vec2D.y = mesh.Coordinate2D[elementID - 1].y[(LocalEdgeNO + 1) % 3] - mesh.Coordinate2D[elementID - 1].y[LocalEdgeNO];
            T norm_e = sqrt(Vec2D.x * Vec2D.x + Vec2D.y * Vec2D.y);
            Vec2D.x /= norm_e, Vec2D.y /= norm_e;

            norm_e /= NumParticlesEachEdge[i] + 1;

            for (uint j = 1; j <= NumParticlesEachEdge[i]; ++j)
            {
                tmpPVEC[j - 1].ElementID = elementID;
                cuDFNsys::Vector2<T> position_rr;
                position_rr.x = mesh.Coordinate2D[elementID - 1].x[LocalEdgeNO] + Vec2D.x * norm_e * j,
                position_rr.y = mesh.Coordinate2D[elementID - 1].y[LocalEdgeNO] + Vec2D.y * norm_e * j;
                tmpPVEC[j - 1].Position2D = position_rr;
            }

            thrust::copy(tmpPVEC.begin(),
                         tmpPVEC.end(),
                         this->ParticlePlumes.begin() + TmpCountParticle);
            TmpCountParticle += NumParticlesEachEdge[i];
        }
    };

    void IdentifyNeighborElements(cuDFNsys::Mesh<T> mesh)
    {
        this->NeighborEleOfOneEle.resize(mesh.Element3D.size());

        std::unordered_map<uint, cuDFNsys::NeighborEle> SharedNodeEle;
        SharedNodeEle.reserve(mesh.Coordinate3D.size());

        //cout << "t 1\n";
        for (uint i = 0; i < mesh.Element3D.size(); ++i)
        {
            uint *Ptr = &(mesh.Element3D[i].x);

            for (uint j = 0; j < 3; ++j)
            {
                //cout << "Ptr[j]: " << Ptr[j] << ", mesh.Coordinate3D.size() = " << mesh.Coordinate3D.size() << endl;
                uint nodeID = Ptr[j] - 1;
                //cout << "\t\t~1" << endl;
                SharedNodeEle[nodeID].NumNeighborEle++;
                if (SharedNodeEle[nodeID].NumNeighborEle > _NumOfNeighborEleAtMost)
                {
                    string AS = "The number of elements sharing the same node exceeds the allowed value\n";
                    throw cuDFNsys::ExceptionsPause(AS);
                }
                //cout << "\t\t~2" << endl;
                SharedNodeEle[nodeID].EleID[SharedNodeEle[nodeID].NumNeighborEle - 1] = i + 1;
                //cout << "\t\t~3" << endl;
                //cout << "\t" << SharedNodeEle[nodeID].NumNeighborEle << ", " << SharedNodeEle[nodeID].EleID[SharedNodeEle[nodeID].NumNeighborEle - 1] << "\n";
            }
        }

        for (uint i = 0; i < mesh.Element3D.size(); ++i)
        {
            uint *Ptr = &(mesh.Element3D[i].x);

            std::set<uint> KL;

            for (uint j = 0; j < 3; ++j)
            {
                for (uint k = 0; k < SharedNodeEle[Ptr[j] - 1].NumNeighborEle; ++k)
                {
                    if (SharedNodeEle[Ptr[j] - 1].EleID[k] != i + 1)
                        KL.insert(SharedNodeEle[Ptr[j] - 1].EleID[k]);
                }
            }

            this->NeighborEleOfOneEle[i].NumNeighborEle = KL.size();

            if (KL.size() > _NumOfNeighborEleAtMost)
            {
                string AS = "The number of neighbering elements of an element exceeds the allowed value\n";
                throw cuDFNsys::ExceptionsPause(AS);
            }
            //cout << "element ID " << i + 1 << ": ";
            uint oo = 0;
            for (std::set<uint>::iterator j = KL.begin(); j != KL.end(); ++j)
            {
                this->NeighborEleOfOneEle[i].EleID[oo] = *j;
                oo++;
                // cout << *j << (oo == KL.size() ? "\n" : ", ");
            }
        }
    };
};
}; // namespace cuDFNsys