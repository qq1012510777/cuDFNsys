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
#include <cmath>

namespace cuDFNsys
{
class ParticleTransport
{
public:
    int NumParticles = 0;
    thrust::host_vector<cuDFNsys::Particle> ParticlePlumes;
    thrust::host_vector<cuDFNsys::EdgeToEle> EdgesSharedEle;
    thrust::host_vector<cuDFNsys::NeighborEle> NeighborEleOfOneEle;

    uint Dir = 2;

private:
    string ParticlePosition = "ParticlePosition.mat";

public:
    ParticleTransport(unsigned long seed,
                      const int &NumOfParticles,
                      const int &NumTimeStep,
                      float delta_T_,
                      float Dispersion_local,
                      thrust::host_vector<cuDFNsys::Fracture> Fracs,
                      cuDFNsys::Mesh mesh,
                      const cuDFNsys::MHFEM &fem,
                      uint Dir_flow,
                      float outletcoordinate)
    {
        this->Dir = Dir_flow;
        this->IdentifyEdgesSharedEle(mesh);
        this->IdentifyNeighborElements(mesh);
        this->InitilizeParticles(NumOfParticles, mesh, fem);
        this->OutputParticleInfoStepByStep(ParticlePosition, 0,
                                           Fracs, mesh);
        // cout << NumTimeStep << ", " << delta_T_ << ", " << Dispersion_local << " ______ \n";
        this->ParticleMovement(seed, NumTimeStep, delta_T_, Dispersion_local, Fracs, mesh, fem, outletcoordinate);
    };

    void ParticleMovement(unsigned long seed,
                          const int &NumTimeStep,
                          float delta_T_,
                          float Dispersion_local,
                          thrust::host_vector<cuDFNsys::Fracture> Fracs,
                          cuDFNsys::Mesh mesh,
                          const cuDFNsys::MHFEM &fem,
                          float outletcoordinate)
    {
        thrust::device_vector<cuDFNsys::Particle> ParticlePlumes_DEV = this->ParticlePlumes;
        thrust::device_vector<cuDFNsys::Fracture> Fracsvec_DEV = Fracs;
        thrust::device_vector<cuDFNsys::EdgeToEle> EdgesSharedEle_vec_dev = this->EdgesSharedEle;
        thrust::device_vector<uint> ElementFracTag_DEV = mesh.ElementFracTag;
        thrust::device_vector<float> Velocity_sep_edge;
        Velocity_sep_edge.reserve(fem.VelocityNormalScalarSepEdges.rows());
        thrust::device_vector<cuDFNsys::EleCoor> Coordinate2D_Vec_dev = mesh.Coordinate2D;
        thrust::device_vector<cuDFNsys::NeighborEle> NeighborEleOfOneEle_dev = this->NeighborEleOfOneEle;

        Eigen::MatrixXf Vg = fem.VelocityNormalScalarSepEdges.cast<float>();
        float *vc = Vg.data();
        Velocity_sep_edge.insert(Velocity_sep_edge.end(), &vc[0], &vc[fem.VelocityNormalScalarSepEdges.rows()]);

        cuDFNsys::Particle *P_DEV = thrust::raw_pointer_cast(ParticlePlumes_DEV.data());
        cuDFNsys::Fracture *Frac_DEV = thrust::raw_pointer_cast(Fracsvec_DEV.data());
        cuDFNsys::EdgeToEle *EdgesSharedEle_DEV = thrust::raw_pointer_cast(EdgesSharedEle_vec_dev.data());
        uint *EleToFracID_ptr = thrust::raw_pointer_cast(ElementFracTag_DEV.data());
        float *velocity_ptr = thrust::raw_pointer_cast(Velocity_sep_edge.data());
        cuDFNsys::EleCoor *Coordinate2D_Vec_dev_ptr = thrust::raw_pointer_cast(Coordinate2D_Vec_dev.data());
        cuDFNsys::NeighborEle *NeighborEleOfOneEle_dev_ptr = thrust::raw_pointer_cast(NeighborEleOfOneEle_dev.data());

        for (uint i = 1; i <= NumTimeStep; ++i)
        {
            time_t t;
            time(&t);

            cout << "The Step " << i << endl;
            ParticleMovementOneTimeStepGPUKernel<<<this->NumParticles / 256 + 1, 256>>>((unsigned long)t + (unsigned long)(i * 100),
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
                                                                                        mesh.Element3D.size());
            cudaDeviceSynchronize();
            this->ParticlePlumes = ParticlePlumes_DEV;
            this->OutputParticleInfoStepByStep(ParticlePosition, i,
                                               Fracs, mesh);
            // if (ParticlePlumes[0].IfReachOutletPlane == true)
            // {
            //     cout << "step " << i << endl;
            // }
            cout << endl;
        }
    };

    void OutputParticleInfoStepByStep(const string &mat_key,
                                      const uint &StepNO,
                                      thrust::host_vector<cuDFNsys::Fracture> Fracs,
                                      cuDFNsys::Mesh mesh)
    {
        cuDFNsys::MatlabAPI M1;

        string mode = "u";
        if (StepNO == 0)
            mode = "w";

        float *particle_position_3D;
        particle_position_3D = new float[this->NumParticles * 3];
        if (particle_position_3D == NULL)
        {
            string AS = "Alloc error in ParticleTransport::OutputParticleInfoStepByStep\n";
            throw cuDFNsys::ExceptionsPause(AS);
        }

        for (size_t i = 0; i < this->NumParticles; ++i)
        {
            float3 tmpPos = make_float3(this->ParticlePlumes[i].Position2D.x,
                                        this->ParticlePlumes[i].Position2D.y,
                                        0);

            uint elementID = this->ParticlePlumes[i].ElementID; // from 1
            uint FracID = mesh.ElementFracTag[elementID - 1];   // from 0

            float Rotate2DTo3D[3][3];
            Fracs[FracID].RoationMatrix(Rotate2DTo3D, 23);

            tmpPos = cuDFNsys::ProductSquare3Float3(Rotate2DTo3D, tmpPos);

            tmpPos = make_float3(tmpPos.x + Fracs[FracID].Center.x,
                                 tmpPos.y + Fracs[FracID].Center.y,
                                 tmpPos.z + Fracs[FracID].Center.z);

            particle_position_3D[i] = tmpPos.x;
            particle_position_3D[i + this->NumParticles] = tmpPos.y;
            particle_position_3D[i + this->NumParticles * 2] = tmpPos.z;
        }
        M1.WriteMat(mat_key, mode, this->NumParticles * 3,
                    this->NumParticles, 3, particle_position_3D, "particle_position_3D_step_" + to_string(StepNO));
        delete[] particle_position_3D;
        particle_position_3D = NULL;

        uint Step[1] = {StepNO};
        M1.WriteMat(mat_key, "u", 1,
                    1, 1, Step, "NumOfSteps");
    };

    void MatlabPlot(const string &mat_key,
                    const string &command_key,
                    const cuDFNsys::Mesh &mesh,
                    const cuDFNsys::MHFEM &fem,
                    const float &L)
    {
        cuDFNsys::MatlabAPI M1;

        size_t node_num = mesh.Coordinate3D.size();

        float *ptr_coordinates_3D;
        ptr_coordinates_3D = new float[node_num * 3];

        if (ptr_coordinates_3D == NULL)
        {
            string AS = "Alloc error in ParticleTransport::MatlabPlot\n";
            throw cuDFNsys::ExceptionsPause(AS);
        }

        for (size_t i = 0; i < node_num; ++i)
        {
            ptr_coordinates_3D[i] = mesh.Coordinate3D[i].x;
            ptr_coordinates_3D[i + node_num] = mesh.Coordinate3D[i].y;
            ptr_coordinates_3D[i + node_num * 2] = mesh.Coordinate3D[i].z;
        }
        M1.WriteMat(mat_key, "w", node_num * 3,
                    node_num, 3, ptr_coordinates_3D, "coordinate_3D");

        delete[] ptr_coordinates_3D;
        ptr_coordinates_3D = NULL;

        //---------------------
        size_t ele_num = mesh.Element3D.size();
        float *ptr_element_3D = new float[ele_num * 3];
        if (ptr_element_3D == NULL)
        {
            string AS = "Alloc error in ParticleTransport::MatlabPlot\n";
            throw cuDFNsys::ExceptionsPause(AS);
        }

        for (size_t i = 0; i < ele_num; ++i)
        {
            ptr_element_3D[i] = mesh.Element3D[i].x;
            ptr_element_3D[i + ele_num] = mesh.Element3D[i].y;
            ptr_element_3D[i + ele_num * 2] = mesh.Element3D[i].z;
        }
        M1.WriteMat(mat_key, "u", ele_num * 3,
                    ele_num, 3, ptr_element_3D, "element_3D");

        delete[] ptr_element_3D;
        ptr_element_3D = NULL;
        //--------------------------

        float *pressure_ELEs = new float[ele_num];
        if (pressure_ELEs == NULL)
        {
            string AS = "Alloc error in ParticleTransport::MatlabPlot\n";
            throw cuDFNsys::ExceptionsPause(AS);
        }

        std::copy(fem.PressureEles.data(),
                  fem.PressureEles.data() + ele_num,
                  pressure_ELEs);

        M1.WriteMat(mat_key, "u", ele_num,
                    ele_num, 1, pressure_ELEs,
                    "pressure_eles");

        delete[] pressure_ELEs;
        pressure_ELEs = NULL;
        //----------------------------------

        std::ofstream oss(command_key, ios::out);
        oss << "clc;\nclose all;\nclear all;\n";
        oss << "load('" << mat_key << "');\n";
        oss << "L = 0.5 * " << L << ";\n";
        oss << "cube_frame = [-L, -L, L; -L, L, L; L, L, L; L -L, L; -L, -L, -L; -L, L, -L; L, L, -L; L -L, -L; -L, L, L; -L, L, -L; -L, -L, -L; -L, -L, L; L, L, L; L, L, -L; L, -L, -L; L, -L, L; L, -L, L; L, -L, -L; -L, -L, -L; -L, -L, L; L, L, L; L, L,-L; -L, L, -L; -L,L, L];\n";
        oss << "figure(1); view(3); title('DFN flow (mhfem)'); xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); hold on\n";
        oss << "patch('Vertices', cube_frame, 'Faces', [1, 2, 3, 4;5 6 7 8;9 10 11 12; 13 14 15 16], 'FaceVertexCData', zeros(size(cube_frame, 1), 1), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0); hold on\n";
        oss << endl;
        oss << "patch('Vertices', coordinate_3D, 'Faces', element_3D, 'FaceVertexCData', pressure_eles, 'FaceColor', 'flat', 'EdgeAlpha', 0.9, 'facealpha', 1); colorbar; view(3); hold on\n";
        oss << "caxis([" << fem.OutletP << ", " << fem.InletP << "]);\n";
        oss << "xlim([-(0.1 * L + L), (0.1 * L + L)])\nylim([ -(0.1 * L + L), (0.1 * L + L) ])\nzlim([ -(0.1 * L + L), (0.1 * L + L) ]);hold on\n\n";
        oss << "S = load('" << ParticlePosition << "');\n";
        oss << "Matrx3D_pso = [];\n";
        oss << "for i = 0:S.NumOfSteps\n";
        //oss << "\ttitle(['DFN flow (mhfem); step NO = ', num2str(i)]);\n";
        oss << "\teval(['Matrx3D_pso(:, :, i + 1) = S.particle_position_3D_step_', num2str(i), '(:, :);']);\n";
        //oss << "\teval(['p_s = scatter3(S.particle_position_3D_step_', num2str(i), '(:, 1), S.particle_position_3D_step_', num2str(i), '(:, 2), S.particle_position_3D_step_', num2str(i), '(:, 3), ''k'', ''o'', ''filled'');']);\n";
        //oss << "\tif (i == 0); pause; else; pause(0.01); end;\n";

        //oss << "\tif (i ~= S.NumOfSteps); delete(p_s); end\n";
        oss << "end\n";
        oss << "N_ = S.NumOfSteps; clear S\n\n";
        oss << "for i = 0:N_\n";
        oss << "\ttitle(['DFN flow (mhfem); step NO = ', num2str(i)]);\n";
        oss << "\tp_s = scatter3(Matrx3D_pso(:, 1, i + 1), Matrx3D_pso(:, 2, i + 1), Matrx3D_pso(:, 3, i + 1), 'k', 'o', 'filled');\n";
        oss << "\tif (i == 0); pause; else; pause(0.01); end;\n";
        oss << "\tif (i ~= N_); delete(p_s); end\n";
        oss << "end\n";

    };

private:
    void IdentifyEdgesSharedEle(cuDFNsys::Mesh mesh)
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
                            cuDFNsys::Mesh mesh,
                            const cuDFNsys::MHFEM &fem){
        //this->ParticlePlumes.resize(NumOfParticles);
        //////////////
        // this->ParticlePlumes.resize(1);
        // this->NumParticles = 1; //mesh.InletEdgeNOLen.size();
        //
        // for (uint i = 8; i < 9; ++i)
        // {
        //     uint EdgeNO = (uint)mesh.InletEdgeNOLen[i].x;     // from 1
        //     uint elementID = (EdgeNO - 1) / 3 + 1;            // from 1
        //     uint FracID = mesh.ElementFracTag[elementID - 1]; // from 0
        //     uint LocalEdgeNO = (EdgeNO - 1) % 3;
        //
        //     float2 position_ = make_float2(0.5f * (mesh.Coordinate2D[elementID - 1].x[LocalEdgeNO] + mesh.Coordinate2D[elementID - 1].x[(LocalEdgeNO + 1) % 3]),
        //                                    0.5f * (mesh.Coordinate2D[elementID - 1].y[LocalEdgeNO] + mesh.Coordinate2D[elementID - 1].y[(LocalEdgeNO + 1) % 3]));
        //
        //     cuDFNsys::Particle tmpP;
        //     tmpP.ElementID = elementID;
        //     tmpP.Position2D = position_;
        //     ParticlePlumes[0] = tmpP;
        // };
        // return;
        ////

        thrust::host_vector<uint> NumParticlesEachEdge(mesh.InletEdgeNOLen.size());
        
        for (size_t i = 0; i < mesh.InletEdgeNOLen.size(); ++i)
        {
            uint EdgeNO = (uint)mesh.InletEdgeNOLen[i].x;
            float length_ = mesh.InletEdgeNOLen[i].y;
        
            float proportion_ = abs(fem.VelocityNormalScalarSepEdges(EdgeNO - 1, 0)) * length_ / fem.QIn;
        
            NumParticlesEachEdge[i] = (uint)round(proportion_ * NumOfParticles);
        }
        
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
        
            float2 position_ = make_float2(0.5f * (mesh.Coordinate2D[elementID - 1].x[LocalEdgeNO] + mesh.Coordinate2D[elementID - 1].x[(LocalEdgeNO + 1) % 3]),
                                           0.5f * (mesh.Coordinate2D[elementID - 1].y[LocalEdgeNO] + mesh.Coordinate2D[elementID - 1].y[(LocalEdgeNO + 1) % 3]));
        
            cuDFNsys::Particle tmpP;
            tmpP.ElementID = elementID;
            tmpP.Position2D = position_;
        
            thrust::host_vector<cuDFNsys::Particle> tmpPVEC(NumParticlesEachEdge[i], tmpP);
        
            thrust::copy(tmpPVEC.begin(),
                         tmpPVEC.end(),
                         this->ParticlePlumes.begin() + TmpCountParticle);
            TmpCountParticle += NumParticlesEachEdge[i];
        }
    };

    void IdentifyNeighborElements(cuDFNsys::Mesh mesh)
    {
        this->NeighborEleOfOneEle.resize(mesh.Element3D.size());

        std::unordered_map<uint, cuDFNsys::NeighborEle> SharedNodeEle;
        SharedNodeEle.reserve(mesh.Coordinate3D.size());

        //cout << "1\n";
        for (uint i = 0; i < mesh.Element3D.size(); ++i)
        {
            uint *Ptr = &(mesh.Element3D[i].x);

            for (uint j = 0; j < 3; ++j)
            {
                //cout << "Ptr[j]: " << Ptr[j] << ";\n";
                SharedNodeEle[Ptr[j] - 1].NumNeighborEle++;
                SharedNodeEle[Ptr[j] - 1].EleID[SharedNodeEle[Ptr[j] - 1].NumNeighborEle - 1] = i + 1;
                //cout << "\t" << SharedNodeEle[Ptr[j] - 1].NumNeighborEle << "\n";
            }
        }
        //cout << "2\n";
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